import psycopg2
import platform
import os.path
import sys
from DatabaseHelper import DatabaseHelper
from MetEngDatabase import MetEngDatabase
from MetEngDatabase import Paper

NOT_PRESENT_NUMBER = -999.0
NOT_PRESENT_STRING = 'NA-Not Avaliable'

DATABASE_DIR   = os.getcwd() + os.sep + 'database' + os.sep
STATIC_DIR     = os.getcwd() + os.sep + 'static' + os.sep
inputTerms = 'InputFields.txt'

inputs = open(STATIC_DIR + inputTerms,"r")

dataTerms = []
muttTerms = []
geneTerms = []

#basically things that are associated with each paper, but not displayed anywhere.
hiddenTerms = []

termDisplayMap  = {}
termCategoryMap = {}
termTypeMap     = {}
termRequiredMap = {}
termEntryMap    = {}
termMultipleMap = {}
parserRules     = {}

for line in inputs:

        line = line.translate(None,'\n\"')
        tokens = line.split("\t")

        category  = tokens[0]
        entryName = tokens[1]
        display   = tokens[2]
        inputType = tokens[3]
        inputRequired = tokens[4]
        multipleValues= tokens[5]
        inputEntry    = tokens[6]
        parsingRule   = tokens[7]

        if(category == 'paper'):
                dataTerms.append(entryName)

        if(category == 'mutant'):
                muttTerms.append(entryName)

        if(category == 'gene'):
                geneTerms.append(entryName)

        if(category == 'hidden'):
                dataTerms.append(entryName)
                hiddenTerms.append(entryName)
        
        termCategoryMap[entryName] = display
        termTypeMap[entryName]     = inputType
        termRequiredMap[entryName] = inputRequired
        termEntryMap[entryName]    = inputEntry
        termMultipleMap[entryName] = multipleValues 
        parserRules[entryName]     = parsingRule

records = os.listdir(DATABASE_DIR)
numMutantsKey = 'NumberofMutants'
numMutationsKey = 'NumberofMutations'

###want to define postgres types based on input file.
#type map from input file
TEXT = 0
NUMBER = 1
DROPDOWN = 2

SINGLE_VALUE = 'single split'
MULTIPLE_VALUE = 'split split'

##postgresTypeMap = {(TEXT, SINGLE_VALUE): 'text',
##                   (TEXT, MULTIPLE_VALUE): 'text[]',
##                   (NUMBER, SINGLE_VALUE): 'real',
##                   (NUMBER, MULTIPLE_VALUE): 'real[]',
##                   (DROPDOWN, SINGLE_VALUE): 'text',
##                   (DROPDOWN, MULTIPLE_VALUE): 'text[]'}

postgresTypeMap = {(TEXT, SINGLE_VALUE): 'text',
                   (TEXT, MULTIPLE_VALUE): 'text',
                   (NUMBER, SINGLE_VALUE): 'text',
                   (NUMBER, MULTIPLE_VALUE): 'text',
                   (DROPDOWN, SINGLE_VALUE): 'text',
                   (DROPDOWN, MULTIPLE_VALUE): 'text'}


paperTable = '(id serial PRIMARY KEY'

#construct create table command.
for term in dataTerms:

    variableType = int(termTypeMap[term])
    valuesType   = parserRules[term]

    term = term.replace(" ","")

    postgresTuple = (variableType, valuesType)

    paperTable = paperTable + ', ' + term + ' ' + postgresTypeMap[postgresTuple]

paperTable = paperTable + ')'

mutantTable = '(id serial PRIMARY KEY,DOI text,mutantNumber integer'

for term in muttTerms:

    variableType = int(termTypeMap[term])
    valuesType   = parserRules[term]

    postgresTuple = (variableType, valuesType)

    mutantTable = mutantTable + ', ' + term + ' ' + postgresTypeMap[postgresTuple]

mutantTable = mutantTable + ')'

geneTable = '(id serial PRIMARY KEY,DOI text,mutantNumber integer, mutationNumber integer'

for term in geneTerms:

    variableType = int(termTypeMap[term])
    valuesType   = parserRules[term]

    postgresTuple = (variableType, valuesType)

    geneTable = geneTable + ', ' + term + ' ' + postgresTypeMap[postgresTuple]

geneTable = geneTable + ')'

annoTable = '(id serial PRIMARY KEY,DOI text,mutantNumber integer, mutationNumber integer, mutation text, annotation text)'

parserRules['mutantNumber'] = 'single split'
parserRules['mutationNumber'] = 'single split'
termTypeMap['mutantNumber'] = '1'
termTypeMap['mutationNumber'] = '1'

try:
        if(platform.system() != 'Windows'):
                #server
                connect = psycopg2.connect("dbname='laser' user='jawi3277' host='localhost' password='winkler'")
        else:
                #client
                connect = psycopg2.connect("dbname='laser' user='james' host='localhost' password='winkler'")
except:
    print "I am unable to connect to the database"
    sys.exit(0)

#controls database access?
cur = connect.cursor()
cur.execute("DROP TABLE IF EXISTS paper")
cur.execute("DROP TABLE IF EXISTS mutant")
cur.execute("DROP TABLE IF EXISTS gene")
cur.execute("DROP TABLE IF EXISTS annotation")

cur.execute("CREATE TABLE paper " + paperTable + ";")
cur.execute("CREATE TABLE mutant " + mutantTable + ";")
cur.execute("CREATE TABLE gene " + geneTable + ";")
cur.execute("CREATE TABLE annotation " + annoTable + ";")


#upload records into database. make sure to sanatize them first!

database = MetEngDatabase(DATABASE_DIR, records, numMutantsKey, numMutationsKey, dataTerms, muttTerms, geneTerms, parserRules, [])
#split each paper into the appropriate dictionaries for insertion into the database
database.parse()

papers = database.records

helper = DatabaseHelper(NOT_PRESENT_NUMBER, NOT_PRESENT_STRING)

for paper in papers:

        paperBacking = paper.paperBacking
        mutantBacking= paper.mutantBacking
        geneBacking  = paper.geneBacking
        annotationBacking = paper.annotationBacking

        for term in hiddenTerms:
                paperBacking[term] = paper.backing[term]

        paperSQL, paperVariables = helper.prepareSQLQuery('paper',paperBacking, termTypeMap, parserRules, postgresTypeMap)
        
        cur.execute(paperSQL,paperVariables)

        numMutants = int(paperBacking[numMutantsKey])

        print paperBacking['Title']

        #unpack the data structures for upload into the database.
        for i in range (1, numMutants+1):

                tempMutantMap = {}

                for term in muttTerms:
                        tempMutantMap[term] = mutantBacking[(i,term)]

                #doi is the unique identifier for each paper/mutant/mutation.

                tempMutantMap['mutantNumber'] = i
                tempMutantMap['DOI'] = paperBacking['DOI']

                mutantSQL, mutantVariables = helper.prepareSQLQuery('mutant',tempMutantMap, termTypeMap, parserRules, postgresTypeMap)
                cur.execute(mutantSQL,mutantVariables)

                numMutations = int(tempMutantMap[numMutationsKey])

                for j in range (1, numMutations+1):

                        tempMutationMap = {}

                        for term in geneTerms:
                                tempMutationMap[term] = geneBacking[(i,j,term)]

                        tempMutationMap['mutantNumber'] = i
                        tempMutationMap['mutationNumber'] = j
                        tempMutationMap['DOI'] = paperBacking['DOI']

                        geneSQL, geneVariables = helper.prepareSQLQuery('gene', tempMutationMap, termTypeMap, parserRules, postgresTypeMap)
                        cur.execute(geneSQL,geneVariables)

                        for key in annotationBacking:

                                tempAnnotationMap = {}

                                tempAnnotationMap['mutantNumber'] = i
                                tempAnnotationMap['mutationNumber'] = j
                                tempAnnotationMap['DOI'] = paperBacking['DOI']
                                tempAnnotationMap['mutation'] = key[3]
                                tempAnnotationMap['annotation'] = annotationBacking[key]

                                rules = {}
                                types = {}

                                for tempKey in tempAnnotationMap:
                                        rules[tempKey] = "single split"
                                        types[tempKey] = "0"

                                annotationSQL, annotationVariables = helper.prepareSQLQuery('annotation', tempAnnotationMap, types, rules, postgresTypeMap)
                                cur.execute(annotationSQL,annotationVariables)


        #write record to database
        connect.commit()

connect.close()
