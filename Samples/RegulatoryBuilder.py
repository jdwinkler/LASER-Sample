from collections import defaultdict
from ModelLogger import ModelLogger
from MetEngDatabase import MetEngDatabase
from MetEngDatabase import Paper
import DatabaseUtils
import networkx
import psycopg2
import psycopg2.extras
import copy
import os

INPUT_DIR = os.getcwd() + os.sep + 'inputs'+ os.sep

class RegulatoryBuilder:

    def __init__(self):
        
        #ecoli regulatory networks from regulonDB
        sigma_gene_file = os.getcwd() + os.sep + 'regulondb' + os.sep + 'network_renamedsigma_gene.txt'
        tf_gene_file = os.getcwd() + os.sep + 'regulondb' + os.sep + 'network_tf_gene.txt'
        tf_operon_file = os.getcwd() + os.sep + 'regulondb' + os.sep + 'network_tf_operon.txt'

        #yeast regulatory network from SGD equivalent
        yeast_file = os.getcwd() + os.sep + 'regulondb' + os.sep + 'yeast_regulatory_network.txt'
        #load networks

        ecoli_edges = RegulatoryBuilder.load_regulon_file(sigma_gene_file,'\t')
        ecoli_edges.extend(RegulatoryBuilder.load_regulon_file(tf_gene_file,'\t'))
        ecoli_edges.extend(RegulatoryBuilder.load_operon_file(tf_operon_file,'\t'))
        yeast_edges = RegulatoryBuilder.load_regulon_file(yeast_file,';')

        EcoliG = networkx.DiGraph(ecoli_edges)
        YeastG = networkx.DiGraph(yeast_edges)

        self.mapper = {'Escherichia coli'.upper():EcoliG,'Saccharomyces cerevisiae'.upper():YeastG}

        self.table_order = {}

        self.table_order['Escherichia coli'.upper()] = ['ecgenes','ecproteins']
        self.table_order['Saccharomyces cerevisiae'.upper()] = ['scgenes','scproteins']

        self.OE_mutations = set(['amplified','oe','plasmid','rbs_tuned'])
        self.REP_mutations = set(['rep','antisense'])
        self.DEL_mutations =set(['del','frameshift','is_insertion','indel'])
        self.RANDOM_mutations = set(['aa_snps','nuc_snps'])
        self.ADD_mutations = set(['sensor','regulated'])

        self.OE_bonus = 2
        self.REP_malus = 0.5

    def getBaseNetworks(self):

        return self.mapper['Escherichia coli'.upper()], self.mapper['Saccharomyces cerevisiae'.upper()]

    @staticmethod
    def search(cur, table, field, englishName, bracketName=None):

        basicQuery = None
        variables = tuple()

        if(bracketName != None):
            basicQuery = 'SELECT * from ' + table + ' where UPPER(' + field + ') = UPPER(%s) or UPPER(' + field + ') = UPPER(%s)'
            variables = (englishName, bracketName)
        if(bracketName == None):
            basicQuery = 'SELECT * from ' + table + ' where UPPER(' + field + ') = UPPER(%s)'
            variables = (englishName,)

        cur.execute(basicQuery, variables)
        records = cur.fetchall()

        if(records != None and len(records) > 0):
            return (records, True)
        else:
            return ([],False)

    #regulonDB files are key (regulator) value (regulatee) directionarylity (+/- or both)
    #will reencode yeast files to same format if possible
    @staticmethod
    def load_operon_file(filename, delim):

        import re

        fhandle = open(filename,'r')
        lines = fhandle.readlines()
        fhandle.close()

        edge_array = []

        for line in lines:

            #comment or empty line
            if(line[0] != '#' and len(line.strip()) > 0):
                tokens = line.split(delim)

                regulator = tokens[0].upper().strip()
                targets = tokens[1].upper().strip()
                
                if(len(tokens) > 2):
                    directionality = tokens[2]
                else:
                    directionality = None

                extracted_targets = re.search(r"\[(.*)\]",targets).group(1)

                target_genes = extracted_targets.split(',')

                for gene in target_genes:

                    gene = gene.strip()
                    edge_array.append((regulator, gene, {'weight':1,'directionality':directionality}))

        return edge_array

    @staticmethod
    def load_regulon_file(filename, delim):

        fhandle = open(filename,'r')
        lines = fhandle.readlines()
        fhandle.close()

        counter = 0

        edge_array = []

        for line in lines:

            #comment or empty line
            if(line[0] != '#' and len(line.strip()) > 0):
                tokens = line.split(delim)

                regulator_toks = tokens[0].upper().strip()

                if(',' in regulator_toks):
                    regulators = regulator_toks.split(",")
                else:
                    regulators = [regulator_toks]
                
                target = tokens[1].upper().strip()
                
                #regulonDB files have this, the SC database does not?
                if(len(tokens) > 2):
                    directionality = tokens[2]
                else:
                    directionality = None

                for regulator in regulators:
                    edge_array.append((regulator, target, {'weight':1,'directionality':directionality}))
                    
                counter = counter + 1

        #G.add_edges_from(edge_array)

        return edge_array

    def laser_to_regnet(self, cur, paper, regulatory_map, numMK='NumberofMutants',numGK='NumberofMutations'):

        #load gene names
        #find if in regulatory network
        #modify accordingly
        #if not, check ecgenes/scgenes for corresponding gene name and converted to accepted accession
        #if still haven't found the gene, quit.
        #return altered network + success/failure log

        mutants = int(paper.paperBacking[numMK])


        output_array = []

        unknownMutation = []

        log_files = []
        networks = []

        table_order = self.table_order
        
        #need to avoid checking deleted genes, would be much faster
        for i in range (1, mutants+1):

            log = ModelLogger('mutant-' + str(i))

            foundGenes = []
            missingGenes = []

            seen_already = set()
            
            mutations = int(paper.mutantBacking[(i,numGK)])
            host = paper.mutantBacking[(i,'Species')]

            regnet = None
            if(host.upper() in regulatory_map.keys()):
                regnet = regulatory_map[host.upper()]
                log.add_section('Selected model',[(paper.paperBacking['Title'], i, host)],'Paper, mutant number, selected regulatory network')
            else:
                log.add_section('Model Pairing Failure', [paper.paperBacking['Title']], 'Could not pair mutant with regulatory network')
                continue
            
            for j in range (1, mutations+1):

                name = paper.geneBacking[(i,j,'GeneName')]
                source = paper.geneBacking[(i,j,'GeneSource')]
                mutation_list = paper.geneBacking[(i,j,'GeneMutation')]

                if(name in seen_already):
                    continue
                else:
                    seen_already.add(name)

                if(source.upper() != host.upper()):
                    #only want genes that are native, ignores heterologous regulatory modifications (I suppose)
                    continue

                #if this gene is already in regnet node list, add to foundGenes and quit.
                if(name.upper() in regnet.nodes()):
                    #original name, converted name, list of mutations to implement, tuple location in case additional info is required
                    foundGenes.append((name, name, mutation_list, [i,j]))
                    continue

                #no simple name match, check databases for synonyms I guess...
                if('[' in name):
                    englishName = name[0:name.find('[')]
                    bracketName = name[name.find('[')+1:name.find(']')]
                else:
                    englishName = name
                    bracketName = None

                if(bracketName != None and bracketName.upper() in regnet.nodes()):
                    foundGenes.append((name, bracketName, mutation_list, [i,j]))
                    continue

                if(bracketName != None and englishName.upper() in regnet.nodes()):
                    foundGenes.append((name, englishName, mutation_list, [i,j]))
                    continue

                tables = table_order[source.upper()]

                #try gene synonyms (if avaliable)
                processed_records = []
                for table in tables:

                    (records, found) = RegulatoryBuilder.search(cur, table, 'common_name', englishName, bracketName)

                    #found some records corresponding to this gene name

                    located_synonym = False

                    if('proteins' in table):
                        for record in records:
                            (conv_records, found) = RegulatoryBuilder.search(cur, table.replace('genes','proteins'), 'unique_id', record['gene'], None)
                            processed_records.extend(conv_records)
                    else:
                        processed_records = records
                    
                for record in processed_records:

                    if(record['synonyms'] != None):                        
                        synonyms = record['synonyms'].split(',')
                        located_synonym = False
                        for synonym in synonyms:
                            if(synonym.upper() in regnet.nodes()):
                                foundGenes.append((name, synonym, mutation_list, [i,j]))
                                located_synonym = True
                                break
                        if(located_synonym == True):
                            break

                    if(record['accession_1'] != None):
                        if(record['accession_1'].upper() in regnet.nodes()):
                            foundGenes.append((name, record['accession_1'], mutation_list, [i,j]))
                            located_synonym = True
                            break

                if(located_synonym == False):
                    missingGenes.append((name, host, source))


            log.add_section('LASER genes missing in regulatory network',missingGenes,'Could not find corresponding node')
            log.add_section('LASER genes discovered in regulatory network',foundGenes,'Paired LASER-Regnet Gene')

            
            log.add_data('foundGenes',foundGenes)
            log.add_data('missingGenes',missingGenes)
            log.add_data('species',host)

            #got the network, now manipulate.

            implemented = []

            tempnet = regnet.copy()

            matchedGenes = []

            for gene in foundGenes:

                (laser_name, regnet_name, mutation_list, laser_tuple) = gene

                regnet_name = regnet_name.upper()

                #this protects against annotation errors where genes are entered more than once
                if(regnet_name not in tempnet.nodes()):
                    continue

                matched = False
                
                for mutation in mutation_list:
                    
                    if(mutation in self.DEL_mutations):
                        #delete node (hopefully deletes associated edges...)
                        tempnet.remove_node(regnet_name)

                        implemented.append((laser_name, 'Deleted node'))
                        matched = True
                        break
                    else:

                        #other mutations affect weight of interactions, I guess.
                        #directionality is already stored, but this way we can represent wider changes.
                        #however, I don't know of a good way to represent mutagenize regulators yet.
                        #will examine that in the future.

                        factor = 0

                        if(mutation in self.OE_mutations or mutation in self.REP_mutations):
                            neighbors = networkx.DiGraph.successors(tempnet,regnet_name)
                            
                            if(mutation in self.OE_mutations):
                                factor = self.OE_bonus
                                implemented.append((laser_name, 'Increased all node-node edge weights'))
                            else:
                                factor = self.REP_malus
                                implemented.append((laser_name, 'Decreased all node-node edge weights'))

                            matched = True

                            for n in neighbors:                            
                                tempnet[regnet_name][n]['weight'] = factor * tempnet[regnet_name][n]['weight']

                        if(mutation in self.RANDOM_mutations):
                            implemented.append((laser_name, 'RANDOM (effects)-Currently Unhandled'))

                        if(mutation in self.ADD_mutations):
                            implemented.append((laser_name, 'ADD (new linkage)-Currently unhandled'))

                if(matched):
                    matchedGenes.append((regnet_name, laser_name, mutation_list))
                else:
                    unknownMutation.append((regnet_name, laser_name, mutation_list))

            log.add_section('Mutated regulatory network',implemented,'Changed node/edge properties')
            log.add_data('matchedGenes',matchedGenes)
            log.add_data('unhandledNodes',unknownMutation)
            log.add_data('unhandledMutations',unknownMutation)

            output_array.append((tempnet,log))

        return output_array

    def generateNetwork(self, cur, paper):

        model_log_tuples = self.laser_to_regnet(cur,paper, self.mapper)

        return model_log_tuples


'''
def main():

    try:
        connect = psycopg2.connect("dbname='biocyc' user='james' host='localhost' password='winkler'")
    except:
        print "I am unable to connect to the database"

    cur = connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    regHelper = RegulatoryBuilder()
    papers = DatabaseUtils.getDatabase()

    missing_dict = defaultdict(int)
    counter = 0
    outputList = []
    uniqueFailedPairings = set([])
    header = ('Gene','Missing Count')

    for paper in papers:
        
        logs,networks = regHelper.generateNetwork(cur, paper)
        for log in logs:
            missingGenes = log.get_data('missingGenes')
            for gene in missingGenes:

                (name, host, species) = gene
                missing_dict[name] = missing_dict[name] + 1
                uniqueFailedPairings.add((name, host, species))

        counter = counter + 1

        if(counter % 10 == 0):
            print 'On paper: %i' % counter

    for key in missing_dict:
        outputList.append((key, str(missing_dict[key])))

    DatabaseUtils.writefile(header, outputList, '\t', 'Genes missing from regulatory networks.txt')
    DatabaseUtils.writefile(('Gene Name','Host','Source'), uniqueFailedPairings, '\t', 'Failed pairings.txt')

    connect.close()




main()
'''
    
    
