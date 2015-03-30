#enzymes.col: describes enzyme:reaction relationships
#reactions.dat can also be used in a similar way (probably better to use)
#genes dat: unique IDs for gene products + their synonyms
#proteins dat: unique ids for reactions from gene products
##
##parsing pattern:
##
##Key _ Value (note that value contains the delimiter '_' on many occasions)
##
##parsing rule:
##
##str.split("_",1), where 1 indicates at most 1 split at the first _ occurrence

## record separator: //

from collections import defaultdict
import psycopg2
import platform
import os.path
import sys
import re

def prepareSQLQuery(table, backing, keys, converter):

    fields =  '(' + ', '.join(list(keys)) + ')'
    arguments = []
    variables = []
    
    for key in keys:
        
            arguments.append('(%s)')

            if(key in converter.keys()):
                key = converter[key]

            if(key in backing.keys()):
                variables.append(unicode(backing[key],"UTF-8","replace"))
            else:
                variables.append(None)

    variables = tuple(variables)
    SQL = 'INSERT INTO ' + table + ' '  + fields + ' VALUES (' + ', '.join(arguments) + ')'
    SQL = unicode(SQL, "UTF-8", "replace")

    return SQL, variables

def tableColumns(columnNames):

    output = 'id serial PRIMARY KEY'

    for name in columnNames:
        output = output + ',' + name + ' text'

    return output

def biocycObjectParser(obj, kv_delimiter, general_entry_delimiter, specific_entry_delimiter):

    objectDictionary = {}
    previousKey = ''

    index = 0
    
    for line in obj:

        line = line.translate(None,'\n\r')

        #removes html tags from input database
        line = re.sub('<.*?>', '', line)

        #continuation from previous key ('/' means line break in the file)
        if('/' == line[0]):

            try:
                objectDictionary[previousKey] = objectDictionary[previousKey] + general_entry_delimiter + line.translate(None,'/')
            except:
                for t in obj:
                    print t

                print 'Line: ' + line

                print previousKey
                    
                raw_input("Failure encountered in multiline parser, press enter to continue")

        else:
            # split only at first hyphen
            tokens = line.split(kv_delimiter,1)

            try:
                key = tokens[0].strip().replace("-","_")
                value = tokens[1].strip()
            except:
                print line
                raw_input("Failure encountered in single line parser, press enter to continue")                  
                
            #multiple entries for same key
            if(key in objectDictionary.keys()):
                oldValue = objectDictionary[key]

                if(key in specific_entry_delimiter.keys()):
                    newValue = oldValue + specific_entry_delimiter[key] + value
                else:
                    newValue = oldValue + general_entry_delimiter + value
                objectDictionary[key] = newValue
            else:
                objectDictionary[key] = value

            #reaction database has coefficients stored oddly
            #will result in a list of coefficients matching the react order
            if( (key == 'LEFT' or key == 'RIGHT') ):
                #still more reactants/products to go
                if(index+1 < len(obj) and obj[index+1].find('^COEFFICIENT') == -1):
                    obj.insert(index+1,'^COEFFICIENT - 1')
                #last entry in reaction database
                if(index+1 >= len(obj)):
                    obj.append('^COEFFICIENT - 1')

            previousKey = key

        index = index + 1

    return objectDictionary


def main():

    METACYC_DIR = os.getcwd() + os.sep + 'reaction_metacyc' + os.sep + 'data' + os.sep
    ECOLI_DIR = os.getcwd() + os.sep + 'reaction_ecoli' + os.sep + 'data' + os.sep
    YEAST_DIR = os.getcwd() + os.sep + 'reaction_yeast' + os.sep + 'data' + os.sep
    ECOLIW_DIR = os.getcwd() + os.sep + 'reaction_ecoliw' + os.sep + 'data' + os.sep

    entry_delimiter = ' '
    kv_delimiter = ' - '

    geneDatabase = 'genes.dat'
    proteinDatabase = 'proteins.dat'
    reactionDatabase = 'reactions.dat'
    compoundDatabase = 'compounds.dat'
    enzymeDatabase = 'enzrxns.dat'
    classDatabase = 'classes.dat'

    connect = psycopg2.connect("dbname='biocyc' user='james' host='localhost' password='winkler'")
    cur = connect.cursor()

    files = []

    files.append((METACYC_DIR + reactionDatabase, 'Reactions'))

    
    files.append((METACYC_DIR + geneDatabase, 'Genes'))
    files.append((METACYC_DIR +  proteinDatabase,'Proteins'))
    files.append((METACYC_DIR + compoundDatabase, 'Compounds'))
    files.append((METACYC_DIR + enzymeDatabase, 'Enzymes'))
    files.append((METACYC_DIR + classDatabase, 'Species'))

    files.append((ECOLI_DIR + geneDatabase, 'EcGenes'))
    files.append((ECOLI_DIR +  proteinDatabase,'EcProteins'))
    files.append((ECOLI_DIR + reactionDatabase, 'EcReactions'))
    files.append((ECOLI_DIR + enzymeDatabase, 'EcEnzymes'))

    files.append((YEAST_DIR + geneDatabase, 'ScGenes'))
    files.append((YEAST_DIR +  proteinDatabase,'ScProteins'))
    files.append((YEAST_DIR + reactionDatabase, 'ScReactions'))
    files.append((YEAST_DIR + enzymeDatabase, 'ScEnzymes'))

    files.append((ECOLIW_DIR + geneDatabase, 'EcwGenes'))
    files.append((ECOLIW_DIR +  proteinDatabase,'EcwProteins'))
    files.append((ECOLIW_DIR + reactionDatabase, 'EcwReactions'))
    files.append((ECOLIW_DIR + enzymeDatabase, 'EcwEnzymes'))

    converter = {}
    #name in column names (below) converted to the biocyc equivalent
    converter['COEFFICIENT'] = '^COEFFICIENT'
    converter['RIGHT_HS'] = 'RIGHT'
    converter['LEFT_HS']  = 'LEFT'
    converter['ORPHAN']  = 'ORPHAN?'

    entries = defaultdict(set)
    entries['Genes'] = set(['UNIQUE_ID','ACCESSION_1','TYPES','COMMON_NAME','PRODUCT','CITATIONS','SYNONYMS'])
    entries['Proteins'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','CATALYZES','ABBREV_NAME','COMPONENT_OF','GENE','LOCATIONS','SPECIES','PRODUCT','SYNONYMS'])
    entries['Reactions'] = set(['UNIQUE_ID','COEFFICIENT','REACTION', 'COMMON_NAME','SYNONYMS','EC_NUMBER', 'ENZYMATIC_REACTION', 'ORPHAN','IN_PATHWAY','LEFT_HS','RIGHT_HS','REACTION_DIRECTION','TYPES'])
    entries['Compounds'] = set(['UNIQUE_ID','COMMON_NAME','INCHI','NON_STANDARD_INCHI','SYNONYMS','CHEMICAL_FORMULA','ATOM_CHARGES','REGULATES','HAS_NO_STRUCTURE','COFACTORS_OF'])
    entries['Enzymes'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','ENZYME','KM','REACTION','SYNONYMS','REGULATED_BY','TEMPERATURE_OPT','SPECIFIC_ACTIVITY','PH_OPT','REACTION_DIRECTION'])
    entries['Species'] = set(['UNIQUE_ID','TYPES','COMMON_NAME'])

    entries['EcGenes'] = set(['UNIQUE_ID','ACCESSION_1','TYPES','COMMON_NAME','PRODUCT','CITATIONS','SYNONYMS'])
    entries['EcProteins'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','CATALYZES','ABBREV_NAME','COMPONENT_OF','GENE','LOCATIONS','SPECIES','PRODUCT','SYNONYMS'])
    entries['EcReactions'] = set(['UNIQUE_ID','COEFFICIENT', 'REACTION','COMMON_NAME','SYNONYMS','EC_NUMBER', 'ENZYMATIC_REACTION', 'ORPHAN','IN_PATHWAY','LEFT_HS','RIGHT_HS','REACTION_DIRECTION','TYPES'])
    entries['EcEnzymes'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','ENZYME','KM','REACTION','SYNONYMS','REGULATED_BY','TEMPERATURE_OPT','SPECIFIC_ACTIVITY','PH_OPT','REACTION_DIRECTION'])

    entries['EcwGenes'] = set(['UNIQUE_ID','ACCESSION_1','TYPES','COMMON_NAME','PRODUCT','CITATIONS','SYNONYMS'])
    entries['EcwProteins'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','CATALYZES','ABBREV_NAME','COMPONENT_OF','GENE','LOCATIONS','SPECIES','PRODUCT','SYNONYMS'])
    entries['EcwReactions'] = set(['UNIQUE_ID','COEFFICIENT', 'REACTION', 'COMMON_NAME','SYNONYMS','EC_NUMBER', 'ENZYMATIC_REACTION', 'ORPHAN','IN_PATHWAY','LEFT_HS','RIGHT_HS','REACTION_DIRECTION','TYPES'])
    entries['EcwEnzymes'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','ENZYME','KM','REACTION','SYNONYMS','REGULATED_BY','TEMPERATURE_OPT','SPECIFIC_ACTIVITY','PH_OPT','REACTION_DIRECTION'])

    entries['ScGenes'] = set(['UNIQUE_ID','ACCESSION_1','TYPES','COMMON_NAME','PRODUCT','CITATIONS','SYNONYMS'])
    entries['ScProteins'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','CATALYZES','ABBREV_NAME','COMPONENT_OF','GENE','LOCATIONS','SPECIES','PRODUCT','SYNONYMS'])
    entries['ScReactions'] = set(['UNIQUE_ID', 'COEFFICIENT', 'REACTION','COMMON_NAME','SYNONYMS','EC_NUMBER', 'ENZYMATIC_REACTION', 'ORPHAN','IN_PATHWAY','LEFT_HS','RIGHT_HS','REACTION_DIRECTION','TYPES'])
    entries['ScEnzymes'] = set(['UNIQUE_ID','TYPES','COMMON_NAME','ENZYME','KM','REACTION','SYNONYMS','REGULATED_BY','TEMPERATURE_OPT','SPECIFIC_ACTIVITY','PH_OPT','REACTION_DIRECTION'])

    
    for key in files:

        fileName = key[0]
        fileType = key[1]

        fhandle = open(fileName,'r')
        lines = fhandle.readlines()
        fhandle.close()

        biocycObjects = []

        objectWrapper = []

        for line in lines:

            #ignore comments in biocyc flat files
            if(line[0] != '#' and len(line) > 0):

                #record element, add to objectWrapper array so long as we haven't found the record separator
                if(line[0:2] != "//"):
                    objectWrapper.append(line)

                #reached record delimiter
                else:
                    biocycObjects.append(objectWrapper)
                    objectWrapper = []

        #so at this point, biocycObjects contains a list for every object in the database.
        #want to iterate over it, convert each entry into standard dictionary, then dump it into an SQL database.

        keys = entries[fileType]
        columnDescriptor = tableColumns(keys)
        #iterate over all records

        #clear existing table in SQL database, then create an empty one in its place

        #fix encoding problem for metacyc database.


        print cur.mogrify("DROP TABLE IF EXISTS " + fileType)
        cur.execute("DROP TABLE IF EXISTS " + fileType)
        cur.execute("CREATE TABLE " + fileType + " (" + columnDescriptor + ")")
        
        print 'Started on ' + fileName + ' database'
        for obj in biocycObjects:

            objDict = {}
            if(fileType.find('Genes')>-1 or fileType.find('Proteins')>-1 or fileType.find('Enzymes') > -1 or fileType.find('Species') >-1):
                delimMap = {}
                delimMap['SYNONYMS'] = ','
                objDict = biocycObjectParser(obj, kv_delimiter, entry_delimiter, delimMap)
            if(fileType.find('Reactions') > -1):
                delimMap = {}
                delimMap['^COEFFICIENT'] = ','
                delimMap['LEFT'] = '+'
                delimMap['RIGHT'] = '+'
                delimMap['IN_PATHWAY'] = ','
                delimMap['ENZYMATIC_REACTION'] = ','
                objDict = biocycObjectParser(obj, kv_delimiter, entry_delimiter, delimMap)
                
            if(fileType.find('Compounds') > -1):
                delimMap = {}
                delimMap['CHEMICAL_FORMULA'] = ','
                delimMap['SYNONYMS'] = ','
                delimMap['REGULATES'] = ','
                delimMap['COFACTORS_OF'] = ','
                objDict = biocycObjectParser(obj, kv_delimiter, entry_delimiter, delimMap)

            #okay, need to do some surgery here. I want species ids in all of the protein tables as follows:
            #ecproteins: species = ECOLI
            #ecwproteins: species = ECOLIW
            #proteins: use whatever is there
            #sgproteins: species = YEAST

            if(fileType.lower() == 'ecproteins'):
                objDict['SPECIES'] = 'TAX-511145'
            if(fileType.lower() == 'ecwproteins'):
                objDict['SPECIES'] = 'TAX-566546'
            if(fileType.lower() == 'scproteins'):
                objDict['SPECIES'] = 'TAX-4932'

            SQL, values = prepareSQLQuery(fileType,objDict, keys, converter)
            try:
                cur.execute(SQL,values)
            except:
                print SQL
                print values
                sys.exit(-1)
            
            # general idea would be to insert this data into a postgres database, using only the desired key/value pairs

        print 'Finished parsing provided biocyc database.'

    connect.commit()

    connect.close()

            


    
            

