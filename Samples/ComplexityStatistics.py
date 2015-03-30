import DatabaseUtils
from collections import defaultdict
from MetEngDatabase import MetEngDatabase
from MetEngDatabase import Paper
from TopologyAnalyzer import TopologyAnalyzer
import networkx
import numpy
import pickle
import scipy.stats

#extracts these fields and writes them out to a set file
def mutation_method_effect_table():

    fhandle = open(DatabaseUtils.INPUT_DIR + 'Term Usage.txt','rU')
    lines = fhandle.readlines()
    fhandle.close()

    outputList = []

    methods = []
    effects = []
    mutations = []

    for line in lines[1:]:

        tokens = line.split('\t')

        key = tokens[0]

        if(key == 'Methods'):
            methods.append(tokens[1])
        if(key == 'Effect'):
            effects.append(tokens[1])
        if(key == 'Mutation'):
            mutations.append(tokens[1])

    #write out info

#outputs an adjacency list (essentially) where col 1 = mutation and col2 = associated method
def mutation_method_mapping_table(papers, numMK = 'NumberofMutants', numGK='NumberofMutations'):

    mut_met_map = defaultdict(set)

    for paper in papers:
        mutants = int(paper.paperBacking[numMK])
        for i in range(1, mutants+1):
            #get methods used by this particular mutant
            methods = paper.mutantBacking[(i,'Method')]
            mutations = int(paper.mutantBacking[(i,numGK)])
            for j in range(1, mutations+1):
                mutation_list = paper.geneBacking[(i,j,'GeneMutation')]
                for mutation in mutation_list:
                    for method in methods:
                        mut_met_map[mutation].add(method)
    #output table
    outputList = []

    for key in mut_met_map:
        methods = mut_met_map[key]
        for method in methods:
            outputList.append((key, method))

    DatabaseUtils.writefile(('Mutation','Method'),outputList,'\t',DatabaseUtils.OUTPUT_DIR + 'Mutation-Method Weighting.txt')

def get_distributions(dataset):

    reg_value_dist = defaultdict(list)
    met_value_dist = defaultdict(list)

    uniqueSubsystems = []

    for point in dataset:

        regulatory = point['regulatory']
        metabolic  = point['metabolic']
        year = point['year']

        met_value_dist['density'].append(point['met_clusters'])
        reg_value_dist['density'].append(point['reg_clusters'])

        for node in regulatory:
            
            for method in regulatory[node]:
                reg_value_dist[method].append(regulatory[node][method])

            reg_value_dist['altered'].append(len(regulatory.keys()))

        uniques = set()

        for node in metabolic:
            
            for method in metabolic[node]:
                
                met_value_dist[method].append(metabolic[node][method])

                if(method == 'subsystem'):
                    uniques.add(metabolic[node][method])

            met_value_dist['altered'].append(len(metabolic.keys()))
            
        numUnique = len(uniques)
        uniqueSubsystems.append(numUnique)

    met_value_dist['subsystem'] = uniqueSubsystems

    return met_value_dist, reg_value_dist


def comparepoint(datadict, value_dist):

    pvalues = defaultdict(dict)

    for node in datadict:

        for field in value_dist:

            #we are examining node-based properties here.
            if('altered' in field):
                continue
            elif('subsystem' in field):
                continue
            elif(field in datadict[node]):
                percent = (100.0 - scipy.stats.percentileofscore(value_dist[field], datadict[node][field]))/100.0
                pvalues[node][field] = percent
        
    return pvalues

def comparedist(datadict, halstead, density, value_dist):

    pvalues = {}

    collapsed = defaultdict(list)

    for node in datadict:
        for field in datadict[node]:
            collapsed[field].append(datadict[node][field])

    collapsed['altered'] = len(datadict.keys())
    collapsed['density'] = density
    collapsed['halstead'] = halstead

    percentile_fields = set(['altered','subsystem','density','halstead'])
    distribution_fields = set(['betweenness','eigenvector','degree'])

    for field in collapsed:

        #there is no collection of values here. we want to find the percentile of the number of nodes altered

        if(field in percentile_fields):

            test_value = collapsed[field]

            if('subsystem' in field):
                test_value = len(set(collapsed[field]))

            if(numpy.isnan(test_value)):
                continue

            percent = (100.0 - scipy.stats.percentileofscore(value_dist[field], test_value))/100.0
            pvalues[field] = percent

        elif(field in distribution_fields):

            (D, p) = scipy.stats.ks_2samp(collapsed[field], value_dist[field])
            pvalues[field] = p
        
    return pvalues

def _generate_stats_container(container, year, datapoint):

    for node in datapoint:
        
        metrics = datapoint[node]
        
        data = container[year]
        for field in metrics:
            if(field not in data):
                data[field] = []
            data[field].append(metrics[field])
            
    if('altered' not in container[year]):
        container[year]['altered'] = []
    container[year]['altered'].append(len(datapoint.keys()))

    return container

def _generate_stats_output(container, filename):

    outputList = []

    value_dist = defaultdict(list)

    header = ('Year','Metric','Samples','Mean','Median','Min Value','Max Value','25th Percentile','75th Percentile','90th Percentile','(Avg) Standard Deviation')

    for year in container:

        for field in container[year]:

            #[year][field] is an array of the data for that field

            data_values = container[year][field]

            #skip string analysis
            if(type(data_values[0]) == str):
                continue
            
            mean = numpy.mean(data_values)
            median = numpy.median(data_values)
            samples = len(data_values)
            stdev = numpy.std(data_values)
            min_value = min(data_values)
            max_value = max(data_values)

            #note that the median is the 50th percentile
            percentile25 = numpy.percentile(data_values, 25)
            percentile75 = numpy.percentile(data_values, 75)
            percentile90 = numpy.percentile(data_values, 90)

            outputTuple = (year, field, samples, mean, median, min_value, max_value, percentile25, percentile75, percentile90, stdev)

            outputList.append(outputTuple)

    DatabaseUtils.writefile(header, outputList, '\t', filename)

#generate descriptive statistics for each metric in the dataset
def generate_statistics_dataset(dataset, met_output, reg_output):
    import numpy

    if(met_output == None or reg_output == None):
        return
    
    #key 1: year of paper
    #key 2: property
    #value: array of results.

    #generate container by year
    met_container = defaultdict(dict)
    reg_container = defaultdict(dict)
    for item in dataset:

        year = item['year']
        met_data = item['metabolic']
        reg_data = item['regulatory']

        if('density' not in met_container[year]):
            met_container[year]['density'] = []
        if('density' not in reg_container[year]):
            reg_container[year]['density'] = []
        if('halstead' not in met_container[year]):
            met_container[year]['halstead'] = []
        if('halstead' not in reg_container[year]):
            reg_container[year]['halstead'] = []
        
        met_container = _generate_stats_container(met_container, year, met_data)
        reg_container = _generate_stats_container(reg_container, year, reg_data)

        met_container[year]['density'].append(item['met_clusters'])
        reg_container[year]['density'].append(item['reg_clusters'])

        met_container[year]['halstead'].append(item['halstead'])
        reg_container[year]['halstead'].append(item['halstead'])

    _generate_stats_output(met_container, met_output)
    _generate_stats_output(reg_container, reg_output)

'''
def statistic_ranking(point, met_value_dist, reg_value_dist):

    regulatory = point['regulatory']
    metabolic  = point['metabolic']

    met_pvalues = TopologyAnalyzer.comparedist(metabolic, met_value_dist)
    reg_pvalues = TopologyAnalyzer.comparedist(regulatory, reg_value_dist)

    return met_pvalues, reg_pvalues

'''
        
#dataset: unpickled list generated by process_dataset function
#met_output: file name for topological metrics for
def statistical_ranking_dataset(dataset, met_output, reg_output):

    if(met_output == None or reg_output == None):
        return

    met_value_dist, reg_value_dist = get_distributions(dataset)

    failed = set()

    regulatory_p_output = []
    metabolic_p_output = []

    keys = ('halstead','altered','subsystem','density','degree','eigenvector', 'betweenness')

    header = ('Title','Mutant Number','Host','altered','subsystem','density','degree','eigenvector', 'betweenness')

    for point in dataset:

        regulatory = point['regulatory']
        metabolic  = point['metabolic']

        met_pvalues = comparedist(metabolic, point['halstead'],point['met_clusters'], met_value_dist)
        reg_pvalues = comparedist(regulatory,point['halstead'],point['reg_clusters'], reg_value_dist)
        
        met_tempList = [point['title'], point['mutant'],point['host']]
        reg_tempList = [point['title'], point['mutant'],point['host']]

        for key in keys:
            
            if(key in met_pvalues):
                met_tempList.append(met_pvalues[key])
            else:
                met_tempList.append(99)

            if(key in reg_pvalues):
                reg_tempList.append(reg_pvalues[key])
            else:
                reg_tempList.append(99)

        regulatory_p_output.append(tuple(reg_tempList))
        metabolic_p_output.append(tuple(met_tempList))

    DatabaseUtils.writefile(header, regulatory_p_output,'\t',reg_output)
    DatabaseUtils.writefile(header, metabolic_p_output,'\t',met_output)

#instead of design level complexity, let's analyze the complexity of each individual intervention
#or at least see what the topological metrics look like for each individual intervention.
def node_complexity_SRD(dataset, met_output, reg_output):

    if(met_output == None or reg_output == None):
        return
    failed = set()

    met_value_dist, reg_value_dist = get_distributions(dataset)

    regulatory_p_output = []
    metabolic_p_output = []

    keys = ('degree','eigenvector','betweenness')

    header = ('Host','Node','Count','degree','eigenvector','betweenness')

    seen_metabolic = set()
    seen_regulatory = set()

    met_node_counter = defaultdict(int)
    reg_node_counter = defaultdict(int)

    for point in dataset:

        regulatory = point['regulatory']
        metabolic  = point['metabolic']

        for node in regulatory:
            reg_node_counter[node] = reg_node_counter[node]+1

        for node in metabolic:
            met_node_counter[node] = met_node_counter[node]+1

    for point in dataset:

        regulatory = point['regulatory']
        metabolic  = point['metabolic']

        met_pvalues = comparepoint(metabolic, met_value_dist)
        reg_pvalues = comparepoint(regulatory, reg_value_dist)
        
        for node in met_pvalues:

            met_tempList = [point['host']]

            if(node in seen_metabolic):
                continue
            else:
                seen_metabolic.add(node)
            
            met_tempList.append(node)
            met_tempList.append(met_node_counter[node])
            for key in keys:
                if(key in met_pvalues[node]):
                    met_tempList.append(met_pvalues[node][key])
                else:
                    met_tempList.append(99)

            metabolic_p_output.append(tuple(met_tempList))

        for node in reg_pvalues:

            reg_tempList = [point['host']]

            if(node in seen_regulatory):
                continue
            else:
                seen_regulatory.add(node)

            reg_tempList.append(node)
            reg_tempList.append(node)
            
            for key in keys:
                if(key in reg_pvalues[node]):
                    reg_tempList.append(reg_pvalues[node][key])
                else:
                    reg_tempList.append(99)

            regulatory_p_output.append(tuple(reg_tempList))
        

    DatabaseUtils.writefile(header, regulatory_p_output,'\t',reg_output)
    DatabaseUtils.writefile(header, metabolic_p_output,'\t',met_output)


def network_comparison(networks, dataset, outputFile, networkTarget = 'metabolic'):

    #extract three topological characteristica: betweenness, degree, eigenvector

    dist = defaultdict(dict)

    for host in networks:

        bc = DatabaseUtils.compute_igraph_betweenness(networks[host])
        ec = DatabaseUtils.compute_igraph_evcentrality(networks[host])
        degrees = networkx.degree(networks[host])

        bc_array = []
        ec_array = []

        for node in bc:
            bc_array.append(bc[node])

        for node in ec:
            ec_array.append(ec[node])

        dist[host]['betweenness'] = bc_array
        dist[host]['eigenvector'] = ec_array
        dist[host]['degree'] = ec_array

    fields = ['betweenness', 'eigenvector', 'degree']

    header = ('Title','Host','Mutant Number','BTW','EVC','DEG')

    outputList = []
    for point in dataset:

        data = point[networkTarget]

        collapsed = defaultdict(list)

        for node in data:
            for field in fields:
                collapsed[field].append(data[node][field])

        tempList = [point['title'], point['mutant'],point['host']]

        for field in fields:
            (D, p) = scipy.stats.ks_2samp(collapsed[field], dist[point['host']][field])
            tempList.append(p)

        outputList.append(tuple(tempList))

    DatabaseUtils.writefile(header, outputList, '\t', outputFile)    

#just outputs design statistics, not comparison to others
#design interconnectivity
#number of nodes altered
#product
#yield/titer
#median node degree, betweenness, eigenvector centrality
def output_design_properties(dataset, papers, outputFile, selection = 'metabolic', numGK='NumberofMutations'):

    #re-org papers into DOI-paper object map

    numeric_fields = ['eigenvector','degree','betweenness']
    count_uniques = ['subsystem']
    strict_counts = ['altered']

    if(selection not in set(['metabolic','regulatory'])):
        print 'Invalid data output selection.'
        return

    paper_dict = {}

    for paper in papers:
        doi = paper.paperBacking['DOI']
        paper_dict[doi] = paper

    if(selection == 'metabolic'):
        density = 'met_clusters'
    if(selection == 'regulatory'):
        density = 'reg_clusters'

    header = ('Title',
              'Model Number',
              'Product',
              'Final Titer',
              'Titer Units',
              'Final Yield',
              'Yield Units',
              'Genes Altered',
              'Nodes Altered',
              '# Subsystems Affected',
              'Perceived Complexity',
              'Clusters Affected',
              'Modified Halstead',
              'Median Degree',
              'Median Betweenness',
              'Median EV')

    keys = ['title',
            'mutantNumber',
            'product',
            'titer',
            'titer_units',
            'yield',
            'yield_units',
            'genes_altered',
            'nodes_altered',
            'subsystem',
            'perception',
            'density',
            'halstead',
            'degree',
            'betweenness',
            'eigenvector']

    outputList = []

    for point in dataset:

        combined_data = {}

        mutantNumber = point['mutant']
        paper = paper_dict[point['DOI']]
        perceived_complexity = paper.paperBacking['project_score']
        product = paper.mutantBacking[(mutantNumber, 'TargetMolecule')]
        ptiter =  paper.mutantBacking[(mutantNumber, 'FinalTiter')]
        titer_units = paper.mutantBacking[(mutantNumber, 'TiterUnits')]
        pyield = paper.mutantBacking[(mutantNumber, 'FinalYield')]
        yield_units = paper.mutantBacking[(mutantNumber, 'YieldUnits')]
        mutations = int(paper.mutantBacking[(mutantNumber, numGK)])

        genes = set()

        for i in range(1,mutations+1):
            genes.add(paper.geneBacking[(mutantNumber,i,'GeneName')])

        combined_data['title'] = paper.paperBacking['Title']
        combined_data['mutantNumber'] = mutantNumber
        combined_data['product'] = product
        combined_data['titer'] = ptiter
        combined_data['titer_units'] = titer_units
        combined_data['yield'] = pyield
        combined_data['yield_units'] = yield_units
        combined_data['density'] = point[density]
        combined_data['perception'] = perceived_complexity

        data = point[selection]
        
        for field in numeric_fields:
            values = []
            for node in data:
                if(field in data[node]):
                    values.append(data[node][field])
            combined_data[field] = numpy.median(values)

        for field in count_uniques:
            unique_values = set()
            for node in data:
                if(field in data[node]):
                    unique_values.add(data[node][field])
            combined_data[field] = len(unique_values)

        combined_data['nodes_altered'] = len(data.keys())
        combined_data['genes_altered'] = len(genes)
        combined_data['halstead'] = point['halstead']

        tempList = []
        for field in keys:
            tempList.append(combined_data[field])
        outputList.append(tuple(tempList))

    DatabaseUtils.writefile(header, outputList, '\t', outputFile)


#analyze_topology_data(DatabaseUtils.OUTPUT_DIR+'Model and Node Topology-+C R-R Network.pkl')

cur_name = DatabaseUtils.OUTPUT_DIR+'Model and Node Topology-+C R-R Network.pkl'

fhandle = open(cur_name,'rb')
dataset = pickle.load(fhandle)
fhandle.close()

papers = DatabaseUtils.getDatabase()
#mutation_method_mapping_table(papers)

ta = TopologyAnalyzer(calculate_tchars = False)

out = DatabaseUtils.OUTPUT_DIR

output_design_properties(dataset, papers, out+'Metabolic Model Propeties (CURRENCY).txt', selection = 'metabolic')
output_design_properties(dataset, papers, out+'Regulatory Model Propeties (CURRENCY).txt', selection = 'regulatory')

statistical_ranking_dataset(dataset, out+'Metabolic Design Rankings Currency.txt', out+'Regulatory Design Rankings Currency.txt')
generate_statistics_dataset(dataset, out+'Metabolic Topology Metrics Currency.txt',out+'Regulatory Topology Metrics Currency.txt')
node_complexity_SRD(dataset, out+'Metabolic Node P-Value Ranking Currency.txt', out+'Regulatory Node P-Value Ranking Currency.txt')

network_comparison(ta.base_met, dataset, DatabaseUtils.OUTPUT_DIR + 'Metabolic Design-Network Comparison.txt', networkTarget ='metabolic')
network_comparison(ta.base_reg, dataset, DatabaseUtils.OUTPUT_DIR + 'Regulatory Design-Network Comparison.txt', networkTarget ='regulatory')




        
