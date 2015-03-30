from collections import defaultdict
from collections import Counter
from MetEngDatabase import MetEngDatabase
from MetEngDatabase import Paper
from cobra import Model, Reaction, Metabolite
from ModelLogger import ModelLogger
import warnings
import DatabaseUtils
import networkx
import cobra
import traceback
import copy
import psycopg2
import psycopg2.extras
import platform
import os.path
import sys
import re
import time

INPUT_DIR = os.getcwd() + os.sep + 'inputs'+ os.sep

class ReactionBuilder:

    #why don't you have static initializers, python?

    def __init__(self, ecoli_xml = 'iJO1366.xml', yeast_xml = 'yeast5.xml'):

        #cobra stuff

        MODEL_DIR = os.getcwd() + os.sep + "metabolic_models" + os.sep

        if(ecoli_xml is not None):
            self.ecoli_model = cobra.io.read_sbml_model(MODEL_DIR + ecoli_xml)
        else:
            self.ecoli_model = None
        if(yeast_xml is not None):
            self.yeast_model = cobra.io.read_sbml_model(MODEL_DIR + yeast_xml)
        else:
            self.yeast_model = None

        if(ecoli_xml == 'iJO1366.xml'):
            self.ecoli_biomass_rxn = ['Ec_biomass_iJO1366_WT_53p95M','Ec_biomass_iJO1366_core_53p95M']
        else:
            self.ecoli_biomass_rxn = None

        #yeast5 model.
        if(yeast_xml == 'yeast5.xml'):
            self.yeast_biomass_rxn = ['r_1670','r_2110']
        else:
            self.yeast_biomass_rxn = None

        #compound maps associate metacyc compounds with their cobra equivalents
        ecoli_compound_dict = ReactionBuilder.setup_compound_map(INPUT_DIR+'Metacyc-Ecoli Metabolite (Compartments).txt','Escherichia coli')
        yeast_compound_dict = ReactionBuilder.setup_compound_map(INPUT_DIR+'Metacyc-Yeast Metabolite Names Synonym Matching.txt', 'yeast')

        #keys to extract number of mutants, number of mutations per mutant
        self.numMK = 'NumberofMutants'
        self.numGK = 'NumberofMutations'

        #how much to increase the lower bound when encountering an OE/DEL mutation
        self.OE_LOWER_BOUND_BONUS = 100
        self.REP_UPPER_BOUND_MALUS = -100

        #group mutations by effect for implement_all_mutations method, terms taken from term usage file
        self.OE_mutations = set(['amplified','oe','plasmid','rbs_tuned'])
        self.REP_mutations = set(['rep','antisense'])
        self.DEL_mutations =set(['del','frameshift','is_insertion','indel'])

        
        #link host background to model, compound dictionary
        self.model_mapper = {'Escherichia coli'.upper(): (self.ecoli_model, ecoli_compound_dict),
                             'Saccharomyces cerevisiae'.upper(): (self.yeast_model, yeast_compound_dict)}

        #tables: complete list of all database tables to search through
        #order: the order to search through these tables, key is host species, value is array with order
        #paper_fields: the field in a LASER record to search through a table (the key) with
        #database_field: the column name for a given table (usually unique_id or common_name, sometimes EC Number
        #papers: list of paper objects representing the laser database
        #nameDictionary; mapping of species name to metacyc id
        #accessionDictionary:
        (self.tables,
         self.searchOrder,
         self.paper_fields,
         self.database_fields,
         self.name_dict,
         self.accession_dict) = ReactionBuilder.setup_constants()

    @staticmethod
    def setup_constants():

        NAME_DICT = INPUT_DIR+'Species-Metacyc Mapping.txt'
        ACCESS_DICT = INPUT_DIR+'ID-Accession Mapping.txt'

        #contains mapping between gene unique id and b/y number for ecoli and yeast
        accessionDictionary = ReactionBuilder.setup_id_mapping(ACCESS_DICT)

        #species name to metacyc equivalent
        nameDictionary = ReactionBuilder.setup_name_dict(NAME_DICT)

        database_fields = {}
        paper_fields = {}

        #setup order to look through tables
        tables = ['ecgenes','ecproteins','ecwproteins','ecwgenes','scgenes','scproteins','enzymelist','genes','proteins','reactions']
        searchOrder = {}
        searchOrder['Escherichia coli'] = ['ecproteins','ecgenes','ecwgenes','ecwproteins','genes','enzymelist','proteins','reactions']
        searchOrder['Saccharomyces cerevisiae'] = ['scproteins','scgenes','genes','enzymelist','proteins','reactions']
        searchOrder['other'] = ['reactions','enzymelist','proteins','genes']

        for table in tables:
            database_fields[table] = 'common_name'
            paper_fields[table] = 'GeneName'

        paper_fields['reactions'] = 'EcNumber'
        database_fields['reactions'] = 'ec_number'

        return (tables, searchOrder, paper_fields, database_fields, nameDictionary, accessionDictionary)

    
    @staticmethod
    def setup_compound_map(compound_file, species):

        fhandle = open(compound_file, 'r')
        lines = fhandle.readlines()
        fhandle.close()

        compound_dict = {}
        species = species.upper()

        for line in lines[1:]:

            tokens = line.strip().split("\t")

            metacyc_unique_id = tokens[0]
            cobra_unique_name = tokens[2]
            compartment = tokens[4]

            if(species == 'Escherichia coli'.upper()):

                if(compartment == 'c'):
                    compartment = 'cytoplasm'
                if(compartment == 'p'):
                    compartment == 'periplasm'
                if(compartment == 'e'):
                    compartment == 'extracellular'

            if(species == 'yeast'.upper() or species == 'Saccharomyces cerevisiae'.upper()):

                compartment = 'unknown'

            #if there is a match, then couldn't associate metacyc/cobra id
            if("$_None" not in metacyc_unique_id):
                compound_dict[(metacyc_unique_id, compartment)] = cobra_unique_name

        return compound_dict

    @staticmethod
    def setup_id_mapping(mapping_file):

        fhandle = open(mapping_file,'r')
        lines = fhandle.readlines()
        fhandle.close()

        id_mapping = {}

        for line in lines[1:]:

            tokens = line.strip().split("\t")

            unique_id = tokens[0]
            accession = tokens[1]
            id_mapping[unique_id] = accession

        return id_mapping

    @staticmethod
    def setup_name_dict(name_file):

        fhandle = open(name_file,'r')

        lines = fhandle.readlines()

        fhandle.close()

        nameDictionary = defaultdict(dict)

        for line in lines[1:]:

            tokens = line.strip().split("\t")

            table = tokens[0]
            species = tokens[1]
            metacyc_id = tokens[2]

            if(table in nameDictionary):
                tempDict = nameDictionary[table]
                tempDict[species] = metacyc_id
                nameDictionary[table] = tempDict
            else:
                nameDictionary[table] = {species : metacyc_id}
                
        return nameDictionary

    #searches for names that match field in the specified table.
    @staticmethod
    def search(cur, table, field, englishName, bracketName=None):

        basicQuery = None
        variables = tuple()

        if(bracketName != None):
            basicQuery = 'SELECT * from ' + table + ' where UPPER(' + field + ') = UPPER(%s) or UPPER(' + field + ') = UPPER(%s) or UPPER(unique_id) = upper(%s)'
            variables = (englishName, bracketName, englishName)
        elif(bracketName == None):
            basicQuery = 'SELECT * from ' + table + ' where UPPER(' + field + ') = UPPER(%s) or UPPER(unique_id) = UPPER(%s)'
            variables = (englishName,englishName)

        cur.execute(basicQuery, variables)
        
        records = cur.fetchall()

        if(records != None and len(records) > 0):
            return (records, True)
        else:
            return ([],False)

    @staticmethod
    def protein_parser(cur, table, laser_species, protein_equiv, metacyc_tax_id = None, use_genus = True):

        keys = []

        #this is the protein in question, check species
        discovered_tax_id = None

        if (protein_equiv['catalyzes'] != None):
            keys = protein_equiv['catalyzes'].split(" ")
            discovered_tax_id = protein_equiv['species']
        else:
            component_of = protein_equiv['component_of']
            
            (comp_records,found) = ReactionBuilder.search(cur, table.replace("genes","proteins"), 'unique_id', component_of)
            protein_complex = None

            if(len(comp_records) > 0):
                protein_complex = comp_records[0]
                discovered_tax_id = protein_complex['species']
                
            if(protein_complex != None and protein_complex['catalyzes'] != None):
                keys = protein_complex['catalyzes'].split(" ")

        #don't bother checking the species
        if(metacyc_tax_id == None):
            return (keys, False)

        #check species or genus for match. species is exact match, genus we have to pull the tax-class from the species table.
        #only exact matches of species tax id permitted
        if(use_genus == False):
            #will fail if no tax id found
            if(metacyc_tax_id == discovered_tax_id):
                return (keys, True)
            else:
                return (keys, False)
        else:
            #can use species genus, but first check to see if species is an exact match.
            if(metacyc_tax_id == discovered_tax_id):
                return (keys, True)
            elif(discovered_tax_id != None):
                #get class
                cur.execute('select * from species where unique_id = %s', (metacyc_tax_id,))
                user_species = cur.fetchone()

                cur.execute('select * from species where unique_id = %s', (discovered_tax_id,))
                discovered_species = cur.fetchone()

                #same class specified in database, related species
                if(user_species['types'] == discovered_species['types']):
                    return (keys, True)

                #check name equivalence, could be related (although I expect the first check to get that) by name only
                discovered_name = discovered_species['common_name'].upper()

                discovered_tokens = discovered_name.split(" ")

                if(len(discovered_tokens) > 1):
                    truncated_dis_name = discovered_name[0] + ". " + discovered_tokens[1].strip()
                else:
                    truncated_dis_name = discovered_name.upper()

                if(laser_species in truncated_dis_name or laser_species in discovered_name):
                    return (keys, True)

                #this is probably a parsing error on the part of the user, but check the genus names.
                user_genus = user_species['common_name'].split(" ")[0].strip().upper()
                discovered_genus = discovered_species['common_name'].split(" ")[0].strip().upper()

                if(user_genus == discovered_genus):
                    return (keys, True)

        #give up
        return (keys, False)

    @staticmethod
    def metacyc_to_rxn(cur, datapoint, nameDictionary):

        #found genes is a string indicating the source table for the data, record is a dictionary with the column headers
        catalyzedrxns = []
        notFound = []
        found = True
        table = datapoint[0]
        record = datapoint[1]
        
        species = record['LASER_SPECIES']

        metacyc_tax_id = None
        species_associated = False

        if(table.replace("genes","proteins") in nameDictionary):
            tempDict = nameDictionary[table.replace("genes","proteins")]
            if(species.upper() in tempDict):
                metacyc_tax_id = tempDict[species.upper()]
                
        #gene tables: get product ID, get enzrxn id, then reaction id.
        keys = []
        enzyme_table = ''
        field = ''
        if(table.find("genes") > -1):
            product_id = record['product']
            enzyme_table = table.replace("genes","enzymes")

            #multiple product ids
            if(" " in product_id):
                temp_ids = product_id.split(" ")
                product_id = temp_ids[-1]
                
            (pe_records,found) = ReactionBuilder.search(cur, table.replace("genes","proteins"), 'unique_id', product_id)    

            if(len(pe_records) > 0):
                protein_equiv = pe_records[0]
                (keys, species_associated) = ReactionBuilder.protein_parser(cur, table, species, protein_equiv, metacyc_tax_id)
            else:
                notFound.append(datapoint)
                
        if(table.find("proteins") > -1):
            enzyme_table = table.replace("proteins","enzymes")
            (keys, species_associated) = ReactionBuilder.protein_parser(cur, table, species, record, metacyc_tax_id)

            if(keys == []):
                notFound.append(datapoint)

        if(table == "enzymelist"):
            keys = [record['unique_id']]
            enzyme_table = 'enzymes'

        if(table == 'reactions'):
            catalyzedrxns.append(record)

        #this means we can associated the gene with a catalytic protein
        if(keys != []):

            for key in keys:
                
                cur.execute('select * from ' + enzyme_table + ' where upper(unique_id) = upper(%s)', (key,))
                #ids for reactions
                enzrxn_results = cur.fetchall()
                for enzrxn_result in enzrxn_results:

                    reaction = enzrxn_result['reaction']

                    #note: I am eliminating any reactions that are unbalanced (totally generic)
                    
                    cur.execute('select * from reactions where upper(unique_id) = upper(%s)', (reaction,))

                    reaction_set = cur.fetchall()

                    for rxn in reaction_set:

                        #do one check: if you cannot
                        catalyzedrxns.append(rxn)

        return (catalyzedrxns, notFound, species_associated)

    @staticmethod
    def compound_matcher_helper(cur, metabolite_array, location, compoundDict):

        compound_information = {}

        #converted-all metabolites, in original order, converted to cobra format if possible or associated with an unique metacyc id
        converted_metabolites = []

        #non-unique metabolites
        special_metabolites = []

        #missing metabolites (no unique matches in metacyc or cobra metabolites)
        missing_metabolites = []

        not_a_cobra_metabolite = []

        for r in metabolite_array:

            r = r.strip()
            #is a known cobra metabolite
            if((r,location) not in compoundDict):
                cur.execute('select * from compounds where unique_id = %s',(r,))
                not_a_cobra_metabolite.append(r)
                associated_compounds = cur.fetchone()
                if(associated_compounds != None):
                    converted_metabolites.append(r)
                    compound_information[r] = associated_compounds
                else:
                    #handle non-specific metabolites here (like menaquinols).
                    #option 1: generate all possible combinations of reactions.
                    #option 2: punt to the end-user, do nothing but record
                    #option 3: pick the best match in the compound converter provided there is only one.

                    #option 3 and a variation of option 2 sound good
                    number_of_matches = 0
                    matching_key = None

                    missing_metabolites.append(r)

                    #strip off bars.
                    temp_r = r.translate(None, '|').strip().upper()
                    for key in compoundDict:
                        if(temp_r in compoundDict[key].upper() and key[1] == location):
                            number_of_matches = number_of_matches + 1
                            matching_key = key

                    if(number_of_matches == 1):
                        converted_metabolites.append(compoundDict[key])
                        special_metabolites.append(r)
                    else:
                        #failure, cannot finding unique match, probably a generic metabolite.
                        #converted_metabolites.append(r)
                        compound_information[r] = None
                        

            else:

                converted_metabolites.append(compoundDict[(r,location)])

        return (converted_metabolites, not_a_cobra_metabolite, special_metabolites, missing_metabolites, compound_information)

    #compound dict: metacyc id to cobra metabolite name
    @staticmethod
    def convert_to_cobra_compounds(cur, rxn_records, location, compoundDict):

        #move these to a file once testing is complete (end of feb?)

        generic_pairing= {'NAD-P-OR-NOP':'NADH-P-OR-NOP','NADH-P-OR-NOP':'NAD-P-OR-NOP'}

        generic_replacements = {'NAD-P-OR-NOP': ['NAD','NADP'], 'NADH-P-OR-NOP' : ['NADH', 'NADPH']}
        partner_compounds = {'NAD':'NADH','NADP':'NADPH','NADPH':'NADP','NADH':'NAD'}

        #contains metabolite information not found in cobra model
        compound_information = {}

        #no converter, don't process
        if(compoundDict == {}):
            return rxn_records

        outputList = []
        not_cobra_metabolites = []
        mapped_metabolites = []
        missing_metabolites = []
        special_metabolites = []

        for record in rxn_records:

            #check for known generic metabolites, replace them in pairs

            for generic_met in generic_pairing:

                if(generic_met in record['left_hs']):

                    #replace, add additional records if needed.
                    replacements = generic_replacements[generic_met]

                    new_reactants = []
                    new_products  = []

                    for rep_met in replacements:
                        new_reactants.append(record['left_hs'].replace(generic_met,rep_met))
                        new_products.append(record['right_hs'].replace(generic_pairing[generic_met], partner_compounds[rep_met]))

                    record['left_hs'] = new_reactants[0]
                    record['right_hs'] = new_products[0]

                    for i in range(1,len(new_reactants)):

                        record_copy = copy.deepcopy(record)
                        record_copy['left_hs'] = new_reactants[i]
                        record_copy['right_hs'] = new_products[i]
                        rxn_records.append(record_copy)

            reactants = record['left_hs'].split("+")
            products = record['right_hs'].split("+")

            (cobra_reactants, not_a_cobra_metabolite, nonuniques, missing, react_cmpd_info) = ReactionBuilder.compound_matcher_helper(cur, reactants, location, compoundDict)
            mapped_metabolites.extend(cobra_reactants)
            missing_metabolites.extend(missing)
            special_metabolites.extend(nonuniques)
            not_cobra_metabolites.extend(not_a_cobra_metabolite)

            for key in react_cmpd_info:
                compound_information[key] = react_cmpd_info[key]
                
            (cobra_products, not_a_cobra_metabolite, nonuniques, missing, prod_cmpd_info) = ReactionBuilder.compound_matcher_helper(cur, products, location, compoundDict)
            mapped_metabolites.extend(cobra_products)
            missing_metabolites.extend(missing)
            special_metabolites.extend(nonuniques)
            not_cobra_metabolites.extend(not_a_cobra_metabolite)

            for key in prod_cmpd_info:
                compound_information[key] = prod_cmpd_info[key]

            record['left_hs'] = "+".join(cobra_reactants)
            record['right_hs']= "+".join(cobra_products)

            outputList.append(record)

        return (outputList, compound_information, mapped_metabolites, not_cobra_metabolites, special_metabolites, missing_metabolites) 
       
    #converts biocyc gene and protein IDs to their cobra equivalents (hopefully)
    @staticmethod
    def gene_name_helper(table, record, nameConverter, cobra_model):

        #procedure depends on table
        if(table == None):
            return None

        #can get b/y numbers from record if obtained from appropriate gene database
        if(table.find('genes') > -1 and record['unique_id'] in nameConverter):
            #b/y number is stored in ACCESSION_1 column in the record

            if(record['accession_1'] in cobra_model.genes):
                return record['accession_1']
            
        if(table.find('proteins') > -1):
            #this is slightly different, because we need to get the gene id, and then the accession number.

            gene = record['gene']
            if(gene != None and gene in nameConverter):
                #nameconverter maps gene unique ids to their accession_1 numbers

                if(nameConverter[gene] in cobra_model.genes):
                    return nameConverter[gene]

        return None

    @staticmethod
    def add_new_reactions(cobra_model, gene_name, metacyc_rxns, metacyc_compounds):

        id_list = []

        unmatchable_metabolites = []

        for rxn in metacyc_rxns:

            rxn_id = rxn['unique_id']

            id_list.append(rxn_id)

            reactants = rxn['left_hs'].split('+')
            products  = rxn['right_hs'].split('+')

            coefficients = rxn['coefficient'].split(",")
            participants = []

            #for matching coefficients
            participants.extend(reactants)
            participants.extend(products)

            reaction = Reaction(rxn_id)
            reaction.name = "Metacyc_" + rxn_id

            location = 'c'

            #can change reaction bounds if you want
            '''
            
            reaction.upper_bound = 1000.  # This is the default
            reaction.objective_coefficient = 0. # this is the default
            '''
            #create metabolite objects

            metabolite_dict = {}
            coeff_index = 0

            reaction_direct = rxn['reaction_direction']

            #some reactions are not associated with a direction, so set to left-to-right
            if(reaction_direct == None):
                reaction_direct = 'LEFT-TO-RIGHT'

            for comp in participants:
                factor = 1
                if(comp in reactants):
                    factor = -1

                coefficient = 0

                if('n' in coefficients[coeff_index]):
                    #can have a range of metabolites, probably won't be balanced...
                    coefficient = 1
                else:
                    coefficient = int(coefficients[coeff_index])

                #reactants are products, products are reactants
                if('RIGHT-TO-LEFT' in reaction_direct.upper()):
                    factor = factor * -1
                    
                if(comp in metacyc_compounds and metacyc_compounds[comp] != None):
                    #metabolite that does not exist in the model
                    comp_dict = metacyc_compounds[comp]
                    formula = comp_dict['chemical_formula'].translate(None,'() ,')
                    metabolite = Metabolite(comp_dict['unique_id'], formula = formula, name = comp_dict['common_name'], compartment = location)
                    metabolite_dict[metabolite] = coefficient * factor
                elif(comp in cobra_model.metabolites):
                    #exists in model, pull it out with get by id
                    metabolite_dict[cobra_model.metabolites.get_by_id(comp)] = coefficient * factor
                else:
                    #a metabolite that we cannot associate with anything properly. insert dummy values for the user to fix.
                    unmatchable_metabolites.append(comp)
                    metabolite = Metabolite(comp, formula = None, name = comp, compartment = location)
                    metabolite_dict[metabolite] = coefficient * factor
            
                coeff_index = coeff_index + 1

            reaction.add_metabolites(metabolite_dict)

            #lower bound must be < 0 to be reversible
            if('REVERSIBLE' in reaction_direct.upper()):
                reaction.lower_bound = -1000

            #cobra reaction rules have no spaces within a gene name
            altered_name = gene_name.replace(" ","_")

            #i suspect that i will have to add a protein complex field to the entry form.
            reaction.gene_reaction_rule = '(' + altered_name + ')'

            cobra_model.add_reaction(reaction)

        return (cobra_model, id_list, altered_name)

    def implement_all_mutations(self, model_updates, cobra_model, objective, log):

        genes_not_added = []
        genes_processed = []
        reactions_added = []
        gene_name_alterations = {}

        #these are the actual reactions touched by the design.
        matchedGenes = []

        unhandled_nodes = []
        unhandled_mutations =[]

        genesModified = 0
        genesNotModified = 0

        for modification in model_updates:

            #name of gene in question
            gene_name = modification['LASER_NAME']
            
            #name of host
            host_species = modification['LASER_HOST']

            #name of gene source
            gene_species = modification['LASER_SPECIES']

            #host accession number if present and host == source
            gene_accession = modification['COBRA_ACCESSION']

            #list of mutations affecting this gene
            mutation_list = modification['LASER_MUTATIONS']

            #rxns corresponding to this gene (might be [] in the future)
            rxns_to_add = modification['CATALYZED_RXNS']

            #compounds associated with the reactions
            compounds_to_add = modification['COMPOUND_INFORMATION']

            added = False

            #spaces are removed when adding genes to cobra models.
            changed_name = gene_name

            if(rxns_to_add != []):

                if(gene_accession == None or host_species != gene_species):
                    #heterologous gene or a native gene not in the model, same difference
                    (cobra_model, id_list, changed_name) = ReactionBuilder.add_new_reactions(cobra_model, gene_name, rxns_to_add, compounds_to_add)
                    reactions_added.extend(id_list)
                    
                    gene_name_alterations[gene_name] = changed_name

                    if(gene_name != changed_name):
                        log.add_section('Gene names changed for cobra compatibility',[(gene_name, changed_name)],'Changed gene names (original, cobra)')
                else:
                    #this is a gene already in the model. probably shouldn't add it.
                    genes_not_added.append((gene_name, 'already present'))

            #no accession: not already in the cobra model
            #changed name not in: not already added to the model
            #heterologous genes with matched reaction fall under this banner.
            if(gene_accession == None and changed_name not in cobra_model.genes):
                genes_not_added.append((gene_name, 'no accession and not added to model already'))
                #skip this loop iteration, go on to next gene
                genesNotModified = genesNotModified + 1
                continue

            if(gene_accession != None):
                gene_id = gene_accession
            else:
                gene_id = changed_name

            cobra_gene = cobra_model.genes.get_by_id(gene_id)
            genes_processed.append((gene_name, gene_id, mutation_list))

            matched = False

            #no _reaction attribute.
            if(hasattr(cobra_gene, '_reaction')):
                reactions_from_genes = cobra_gene.get_reaction()
            else:
                reactions_from_genes = []

            for mutation in mutation_list:

                #delete gene, break loop.
                if(mutation in self.DEL_mutations):
                    cobra_gene.remove_from_model()
                    matched = True
                    break

                for rxn in reactions_from_genes:
                    if(mutation in self.OE_mutations):
                        #get corresponding reaction
                        rxn.lower_bound = rxn.lower_bound + self.OE_LOWER_BOUND_BONUS
                        matched = True

                    if(mutation in self.REP_mutations):
                        rxn.upper_bound = rxn.upper_bound + self.REP_UPPER_BOUND_MALUS
                        matched = True

            if(matched):
                genesModified = genesModified + 1
                
                for rxn in reactions_from_genes:
                    #print rxn.id, gene_name, gene_id, mutation_list
                    matchedGenes.append((rxn.id, gene_name, mutation_list, gene_id))
            else:
                genesNotModified = genesNotModified + 1
                for rxn in reactions_from_genes:
                    unhandled_nodes.append((rxn.id, gene_name, mutation_list, gene_id))

                #genes that are in the model but have no associated reactions fall here.
                if(reactions_from_genes == []):
                    genes_not_added.append((gene_id, 'No associated reactions'))
                #genes that do not have implementable mutations (at least, currently)
                else:
                    unhandled_mutations.append((gene_id,'Mutation not found'))

        ##debug information
        log.add_section('Genes not added to model',genes_not_added,'Genes not added to model due to')

        log.add_data('not_modified_genes',genesNotModified)
        log.add_data('modified_genes',genesModified)
        log.add_data('matchedGenes',matchedGenes)
        log.add_data('unhandledMutations',unhandled_mutations)
        log.add_data('unhandledNodes',unhandled_nodes)

        failed_data_hook = []
        for gene in genes_not_added:
            if(gene[1] != 'already present'):
                failed_data_hook.append(gene)

        log.add_data('failedGenes',failed_data_hook)
        
        reactions_debug = []

        for item in reactions_added:
            reactions_debug.append((item, cobra_model.reactions.get_by_id(item).reaction))

        log.add_data('reactions_added',reactions_added)
        log.add_data('reactions_debug',reactions_debug)
        
        log.add_section('Final reactions added to cobra model',reactions_debug,'Reaction IDs added to cobra model, actual reaction')
        log.add_section('Processed mutations',genes_processed,'Name, cobra id, and mutations of processed genes')

        #try to get objective reaction

        if(objective['unique_id'] in cobra_model.metabolites):
            metObj = cobra_model.metabolites.get_by_id(objective['unique_id'])

            obj_reactions = []

            for rxn in metObj.get_reaction():
                obj_reactions.append(rxn.reaction)
            log.add_data('foundObjective',True)
            log.add_section('Objective-Reaction Linkage',obj_reactions,'Discovered reaction with %s objective metabolite participating' % objective['unique_id'])
        else:
            log.add_data('foundObjective',False)
            log.add_section('Objective-Reaction Linkage',[],'Discovered reaction with %s objective metabolite participating' % objective)

        #generate_report(cobra_model, paper_title, paper_filename, genes_not_added, genes_processed, reactions_added)

        return cobra_model, log

    def laser_to_cobra(self, cursor, paper, tables, paper_field, database_field, searchOrder, nameDictionary, accessionDictionary, model_mapper, numMK='NumberofMutants',numGK='NumberofMutations'):

        mutants = int(paper.paperBacking[numMK])

        log_files = []
        modified_models = []

        output_array = []

        
        processed_genes = 0
        
        #need to avoid checking deleted genes, would be much faster
        for i in range (1, mutants+1):

            log = ModelLogger('mutant-' + str(i))

            LASER_to_METACYC = 0

            missingGenes = []
            foundGenes = []

            mutations = int(paper.mutantBacking[(i,numGK)])
            host = paper.mutantBacking[(i,'Species')]

            #get reference to cobra model
            model = None
            compoundDict = None
            if(host.upper() in model_mapper.keys()):
                model_tuple = model_mapper[host.upper()]
                model = model_tuple[0]
                compoundDict = model_tuple[1]

                log.add_section('Selected model',[(paper.paperBacking['Title'], i, host)],'Paper, mutant number, selected cobra model')

            #skip this mutant if it has no corresponding model
            if(model == None):
                log.add_section('Model Pairing Failure', [paper.paperBacking['Title']], 'Could not pair mutant with cobra model')
                continue

            check_duplicates_input = set()
            
            for j in range (1, mutations+1):

                if(paper.geneBacking[(i,j,'GeneName')].upper() != 'NONE'):
                    processed_genes = processed_genes + 1

                geneName = paper.geneBacking[(i,j,'GeneName')]            
                mutationList = paper.geneBacking[(i,j,'GeneMutation')]
                source = paper.geneBacking[(i,j,'GeneSource')]
                queryTables = []

                if(geneName in check_duplicates_input):
                    continue
                else:
                    check_duplicates_input.add(geneName)

                if(source in searchOrder):
                    queryTables = searchOrder[source]
                else:
                    queryTables = searchOrder['other']

                finished = False

                #optimizations: don't bother processing genes already in the cobra models-TODO
                #don't process deleted genes-TODO
                
                for table in queryTables:

                    #select table/field to search
                    searchPaperField    = paper_field[table]
                    name = paper.geneBacking[(i,j,searchPaperField)]


                    #skip these entries, usually padding/garbage
                    if(name.upper() == "NONE" or name == ''):
                        continue
                    
                    searchDatabaseField = database_field[table] 
                    if(name.find('[') > -1):
                        englishName = name[0:name.find('[')]
                        bracketName = name[name.find('[')+1:name.find(']')]
                        (results, found) = ReactionBuilder.search(cursor,table,searchDatabaseField,englishName,bracketName)
                    else:
                        englishName = name
                        (results, found) = ReactionBuilder.search(cursor,table,searchDatabaseField,englishName,None)

                    if(found):
                        
                        LASER_to_METACYC = LASER_to_METACYC + 1
                        finished = True
                        
                        for record in results:
                            record['LASER_TUPLE'] = (i,j)
                            record['LASER_MUTATIONS'] = mutationList
                            record['LASER_NAME'] = name
                            record['LASER_HOST'] = host
                            record['LASER_SPECIES'] = source
                            foundGenes.append((table,record))
                        
                        break

                if(not finished):
                    #these are genes that cannot be paired to anything in metacyc
                    missingGenes.append((geneName, source))

            ##at this point, you have laser->metacyc. want metacyc id to rxn.

            #print 'Processed gene count %i, missing + found count % i' % (processed_genes, LASER_to_METACYC + len(missingGenes))

            #don't bother with the missing (including deleted ones) genes when searching for reactions at this stage.

            existingRxns = set()
            existingGenes = set()
            annotatedGenes = []

            speciesAssociated = []
            specifiedLocation = []
            missing_rxns = []
            found_rxns = []

            mapped_metabolites = []
            unpaired_cobra_metabolites = []
            missing_metabolites = []
            special_metabolites = []
            converted_metabolites = []

            check_duplicates_biocyc = set()

            for gene in foundGenes:

                table = gene[0]
                record = gene[1]
                
                #list of reactions catalyzed by this gene
                (initialrxns, missing_gene_rxn, species_associated) = ReactionBuilder.metacyc_to_rxn(cursor, gene, nameDictionary)
                
                #get the accession number of this gene
                accession = ReactionBuilder.gene_name_helper(table, record, accessionDictionary, model)

                #assume cytoplasmic expression unless specified in mutation list

                warnings.warn('Reminder: need to specify correct compartmentalization for all metabolites in yeast, E. coli models.')
                if(host.upper() == 'Escherichia coli'.upper()):
                    location = 'cytoplasm'
                elif(host.upper() == 'Saccharomyces cerevisiae'.upper()):
                    location = 'unknown'
                else:
                    location = 'new organism?'

                if('compartmentalization' in record['LASER_MUTATIONS']):
                    (i,j) = record['LASER_TUPLE'] 
                    location = paper.annotationBacking[(i,j,'GeneMutation','compartmentalization')]


                (catalyzedrxns, compound_information, mapped_mets, non_cobra, nonuniques, missing_mets) = ReactionBuilder.convert_to_cobra_compounds(cursor, initialrxns, location, compoundDict)

                sievedrxns = []
                for rxn in catalyzedrxns:
                    if(rxn['unique_id'] not in existingRxns):
                        sievedrxns.append(rxn)
                        existingRxns.add(rxn['unique_id'])

                record['CATALYZED_RXNS'] = sievedrxns
                record['COBRA_ACCESSION'] = accession
                record['COMPOUND_INFORMATION'] = compound_information

                if(record['LASER_NAME'].upper() not in existingGenes):

                    speciesAssociated.append((record['LASER_NAME'], species_associated))

                    specifiedLocation.append((record['LASER_NAME'], location))

                    #don't care about reactions that get deleted or are already present in the cobra model
                    if('del' not in record['LASER_MUTATIONS'] and accession == None):
                        
                        if(missing_gene_rxn != []):
                            missing_rxns.append(record['LASER_NAME'])

                        for rxn in initialrxns:
                            found_rxns.append((record['LASER_NAME'], rxn['unique_id'], rxn['left_hs'] + "->" + rxn['right_hs']))

                        for met in non_cobra:
                            unpaired_cobra_metabolites.append((record['LASER_NAME'], met))

                        for met in mapped_mets:
                            mapped_metabolites.append((record['LASER_NAME'], met))

                        for met in missing_mets:
                            missing_metabolites.append((record['LASER_NAME'],met))

                        for elem in nonuniques:
                            special_metabolites.append((record['LASER_NAME'], elem))

                        if(compound_information != {}):
                            for key in compound_information:
                                compound = compound_information[key]

                                if(compound != None):
                                    converted_metabolites.append((compound['unique_id'], compound['common_name']))
                                else:
                                    converted_metabolites.append(('No Metacyc ID', key))
                    
                    annotatedGenes.append(record)
                    existingGenes.add(record['LASER_NAME'].upper())
                    

            set_filter = set()
            temp_array = []
            for gene in foundGenes:
                table = gene[0]
                record = gene[1]
                if(record['LASER_NAME'].upper() not in set_filter):
                    temp_array.append((record['LASER_NAME'],record['unique_id'], table))
                    set_filter.add(record['LASER_NAME'].upper())

            #data for database analysis

            log.add_data('missingGenes',missingGenes)
            log.add_data('foundGenes', temp_array)
            log.add_data('title',paper.paperBacking['Title'])

            #report for users

            log.add_section('Initial LASER -> Metacyc Gene Pairing Failures', missingGenes, 'Could not pair LASER to metacyc ID')
            log.add_section('Identified gene-metacyc pairings',temp_array,'LASER Gene Name, unique id in metacyc/biocyc, corresponding table')
            log.add_section('Species-Gene Association',speciesAssociated,'Gene is paired with Metacyc/Biocyc species (Gene, True/False)')
            log.add_section('Missing metacyc ID -> metacyc rxn pair',missing_rxns,'Could not pair metacyc ID to reactions for gene')
            log.add_section('Metacyc ID -> metacyc rxn pairing',found_rxns,'Gene name, unique ID for reaction, formula')
            log.add_section('Reaction Location',specifiedLocation,'Location of metabolites for specified gene/reaction')
            log.add_section('Added metabolites',converted_metabolites,'Metabolites added from reactions')
            log.add_section('Missing metabolites',missing_metabolites,'Metabolites missing from reactions')
            log.add_section('Generic metabolites',special_metabolites,'Generic metabolites in reactions')

            log.add_data('missing_rxns',missing_rxns)
            log.add_data('mapped_metabolites',mapped_metabolites)
            log.add_data('added_metabolites',converted_metabolites)
            log.add_data('missing_metabolites',missing_metabolites)
            log.add_data('generic_metabolites',special_metabolites)
            log.add_data('speciesAssociated',speciesAssociated)
            log.add_data('species',host)
            log.add_data('successful_pairings', LASER_to_METACYC)            

            objective_metabolite = paper.mutantBacking[(i,'TargetMolecule')]

            #try to get compound associated with target molecule
            if(objective_metabolite not in compoundDict):
                cursor.execute('select * from compounds where upper(common_name) = upper(%s) or upper(unique_id) = upper(%s)',(objective_metabolite,objective_metabolite))
                objective_compound = cursor.fetchone()
                if(objective_compound == None):
                    objective_compound = {'unique_id':objective_metabolite}
            else:
                objective_compound = {'unique_id' : compoundDict[objective_metabolite]}

            log.add_section('Target metabolite for production',[objective_metabolite, objective_compound],'Provided and converted name (hopefully)')

            log.add_data('target_metabolite',objective_metabolite)

            #need to add export reaction for objective metabolite at future date

            #matched names is defined in this function
            modified_model, log = self.implement_all_mutations(annotatedGenes, model.copy(), objective_compound, log)


            output_array.append((modified_model,log))

        return output_array

    @staticmethod
    #this is the network structure that Barabisi et al used in the "The large-scale organization of metabolic networks" (Nature, 2000)
    #I personally don't really like it as it creates a massive complex web that is difficult to parse, but eh.
    #reviewers will probably like to see multiple network architectures
    def _met_rxn_bipartite_network(model, ignore_reactions, ignore_metabolites):

        G = networkx.DiGraph()
        reaction_set = model.reactions

        cache = defaultdict(dict)
        already_processed = set()
        skip = set()

        edge_set = []

        weight = {'weight':1.0}

        for bmn in ignore_reactions:
            already_processed.add(bmn)
            skip.add(bmn)

        for rxn in reaction_set:

            isReversible = rxn.reversibility

            #skip these reactions
            if(rxn.id in ignore_reactions):
                continue

            reactants = rxn.reactants
            products = rxn.products

            for reactant in reactants:

                #skips metabolite if excluded from network generation
                if(reactant.id in ignore_metabolites or reactant.name in ignore_metabolites):
                    continue

                #reactant to rxn node
                edge_set.append((reactant.id,rxn.id,weight))
                if(isReversible):
                    #can also walk back from the reactant (i.e. product) to the rxn node
                    edge_set.append((rxn.id, reactant.id, weight))

            for product in products:
                #rxn node to products

                #skips metabolite if excluded from network generation
                if(product.id in ignore_metabolites or product.name in ignore_metabolites):
                    continue
                
                edge_set.append((rxn.id, product.id, weight))
                if(isReversible):
                    #can convert products back (through the reaction node) to reactants
                    edge_set.append((product.id, rxn.id, weight))

        G.add_edges_from(edge_set)
        return G
    
    @staticmethod
    def _rxn_rxn_network(model, ignore_reactions, ignore_metabolites):

        G = networkx.DiGraph()
        reaction_set = model.reactions

        cache = defaultdict(dict)
        already_processed = set()
        skip = set()

        edge_set = []

        weight = {'weight':1.0}

        for bmn in ignore_reactions:
            already_processed.add(bmn)
            skip.add(bmn)
        
        for rxn_source in reaction_set:
            cache[rxn_source.id]['reactants'] = set([rxn.id for rxn in rxn_source.reactants])
            cache[rxn_source.id]['products'] = set([rxn.id for rxn in rxn_source.products])

        #identifies currency metabolites to ignore for rxn-met or rxn-rxn comparison
        exmets = ignore_metabolites

        for rxn_source in reaction_set:

            #want to capture reversibilities in the network
            #isReversible = rxn_source.lower_bound < 0 and rxn_source.upper_bound > 0

            if(rxn_source.id in skip):
                continue

            #inactive reaction, skip
            if(rxn_source.lower_bound == 0 and rxn_source.upper_bound == 0):
                continue

            #get metabolites
            reactants = cache[rxn_source.id]['reactants']
            products = cache[rxn_source.id]['products']

            for rxn_target in reaction_set:

                if(rxn_source.id == rxn_target.id):
                    continue

                if(rxn_target.id in already_processed):
                    continue

                if((rxn_target.lower_bound == 0 and rxn_target.upper_bound == 0)):
                    continue

                #check if reactants == products or products == reactants
                treactants = cache[rxn_target.id]['reactants']
                tproducts = cache[rxn_target.id]['products']
                
                #remove non-exchange metabolite, these sets are any known trivial shared metabolites
                shared_trp = (treactants & products) - exmets
                shared_tpr = (tproducts & reactants) - exmets

                if(len(shared_trp) > 0 and len(treactants) > 0 and len(products) > 0):
                    #print shared_trp
                    edge_set.append((rxn_source.id,rxn_target.id, weight))

                if(len(shared_tpr) > 0 and len(tproducts) > 0 and len(reactants) > 0):
                    edge_set.append((rxn_target.id,rxn_source.id, weight))
                

            already_processed.add(rxn_source.id)

        G.add_edges_from(edge_set)

        return G

    @staticmethod
    #converts a metabolic network in a cobra model to a networkx directed graph
    #nodes: reaction ids, metabolite ids
    #directed (only bidirectional links if reactions is reversible)
    #option to ignore metabolites with high average degrees, like ATP, NAD/NADP, NADH/NADPH, etc
    #option to remove specific reactions by adding them to biomass_node (an iterable of values)
    def extractTopology(model, method, ignore_reactions = set(), ignore_metabolites = set()):

        if(method.upper() == 'RXN-RXN'):
            return ReactionBuilder._rxn_rxn_network(model, ignore_reactions, ignore_metabolites)
        if(method.upper() == 'MET-RXN'):
            return ReactionBuilder._met_rxn_bipartite_network(model, ignore_reactions, ignore_metabolites)

        print 'No valid network generation method specified.'
        return None

    def getBaseModels(self):

        return self.ecoli_model, self.yeast_model
                
    def generateModel(self, cur, paper):

        model_log_tuples = self.laser_to_cobra(cur, paper, self.tables, self.paper_fields, self.database_fields, self.searchOrder, self.name_dict, self.accession_dict, self.model_mapper, numMK = self.numMK, numGK = self.numGK)
        return model_log_tuples

'''
            
def main():

    try:
        connect = psycopg2.connect("dbname='biocyc' user='james' host='localhost' password='winkler'")
    except:
        print "I am unable to connect to the database"

    cur = connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    rxnHelper = ReactionBuilder()

    papers = DatabaseUtils.getDatabase()

    model_tuple = rxnHelper.generateModel(cur, papers[0])

    model_tuple[0][1].print_report()

    print model_tuple[0][1].get_data('matchedGenes')

    connect.close()

main()
'''

