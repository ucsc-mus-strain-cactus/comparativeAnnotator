import os
from itertools import izip_longest, product
import cPickle as pickle

import src.classifiers, src.details, src.attributes
import lib.sqlite_lib as sql_lib
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger


class ConstructDatabases(Target):
    def __init__(self, outDir, genomes, classifiers, details, attributes, alnPslDict, primaryKeyColumn):
        Target.__init__(self)
        self.outDir = outDir
        self.genomes = genomes
        self.classifiers = classifiers
        self.details = details
        self.attributes = attributes
        self.alnPslDict = alnPslDict
        self.primaryKeyColumn = primaryKeyColumn

    def run(self):
        logger.info("Merging pickled files into databases")
        classifyDb = os.path.join(self.outDir, "classify.db")
        detailsDb = os.path.join(self.outDir, "details.db")
        attributesDb = os.path.join(self.outDir, "attributes.db")
        self.initializeDb(classifyDb)
        self.initializeDb(detailsDb)
        self.initializeDb(attributesDb)
        for classifier, genome in product(self.classifiers, self.genomes):
            valueDict = pickle.load(open(self.globalTempDir(), classifier.__name__ + genome, "rb"))
            self.simpleUpdateWrapper(valueDict)
        for detail, genome in product(self.details, self.genomes):
            valueDict = pickle.load(open(self.globalTempDir(), detail.__name__ + genome, "rb"))
            self.simpleBedUpdateWrapper(valueDict)
        for attribute, genome in product(self.attributes, self.genomes):
            valueDict = pickle.load(open(self.globalTempDir(), attribute.__name__ + genome, "rb"))
            self.simpleUpdateWrapper(valueDict)

    def invertDict(self, d):
        for a, b in d.iteritems():
            yield b, a

    def simpleUpdateWrapper(self, valueDict):
        """
        If your classifier is going to do a simple 1-1 update with a valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            sql_lib.updateRows(cur, self.genome, self.primaryKey, self.getColumn(), 
                    self.invertDict(valueDict))

    def simpleBedUpdateWrapper(self, valueDict):
        """
        If your details-mode classifier has BED records for values in its valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            sql_lib.updateRows(cur, self.genome, self.primaryKey, self.getColumn(), 
                    self.detailsEntryIter(valueDict.iteritems()))     

    def detailsEntryIter(self, valueIter):
        """
        General case for converting a valueIter to a string entry representing 1 or more BED records
        valueIters are generally lists of lists where each sublist represents a BED record
        """
        for aId, entry in valueIter:
            if entry is None:
                yield None, aId
            elif type(entry[0]) != list:
                #only one entry
                yield "\t".join(map(str,entry)), aId
            else:
                bedEntries = ["\t".join(map(str, x)) for x in entry]
                yield "\n".join(bedEntries), aId

    def initializeDb(self, dbPath):
        columnDefinitions = [[x.__name__, x._getType()] for x in self.classifiers]
        #find alignment IDs from PSLs (primary key for database)
        for genome in self.genomes:
            aIds = set(x.split()[9] for x in open(self.alnPslDict[aId]))
            initializeSqlTable(dbPath, genome, columnDefinitions, self.primaryKeyColumn)
            initializeSqlRows(dbPath, genome, aIds, self.primaryKeyColumn)

    def initializeSqlTable(self, db, genome, columns, primaryKey):
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.initializeTable(cur, genome, columns, primaryKey)

    def initializeSqlRows(self, db, genome, aIds, primaryKey):
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.insertRows(cur, genome, primaryKey, [primaryKey], izip_longest(aIds, [None]))