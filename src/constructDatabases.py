import os
from itertools import izip_longest, product
import cPickle as pickle

import src.classifiers, src.attributes
import lib.sqlite_lib as sql_lib
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger


class ConstructDatabases(Target):
    def __init__(self, outDir, genomes, classifiers, attributes, alnPslDict, primaryKeyColumn):
        Target.__init__(self)
        self.outDir = outDir
        self.genomes = genomes
        self.classifiers = classifiers
        self.attributes = attributes
        self.alnPslDict = alnPslDict
        self.primaryKeyColumn = primaryKeyColumn

    def run(self):
        logger.info("Merging pickled files into databases")
        classifyDb = os.path.join(self.outDir, "classify.db")
        detailsDb = os.path.join(self.outDir, "details.db")
        attributesDb = os.path.join(self.outDir, "attributes.db")
        self.initializeDb(classifyDb, self.classifiers, dataType="INTEGER")
        self.initializeDb(detailsDb, self.classifiers, dataType="TEXT")
        for classifier, genome in product(self.classifiers, self.genomes):
            classifyDict = pickle.load(open(os.path.join(self.outDir, genome, "Classify" + classifier.__name__ + genome), "rb"))
            self.simpleUpdateWrapper(classifyDict, classifyDb, genome, classifier.__name__)
            detailsDict = pickle.load(open(os.path.join(self.outDir, genome, "Details" + classifier.__name__ + genome), "rb"))
            self.simpleBedUpdateWrapper(detailsDict, detailsDb, genome, classifier.__name__)
        self.initializeDb(attributesDb, self.attributes)
        for attribute, genome in product(self.attributes, self.genomes):
            attributeDict = pickle.load(open(os.path.join(self.outDir, genome, "Attribute" + attribute.__name__ + genome), "rb"))
            self.simpleUpdateWrapper(attributeDict, attributesDb, genome, attribute.__name__)

    def invertDict(self, d):
        for a, b in d.iteritems():
            yield b, a

    def simpleUpdateWrapper(self, valueDict, db, genome, column):
        """
        If your classifier is going to do a simple 1-1 update with a valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.updateRows(cur, genome, self.primaryKeyColumn, column, self.invertDict(valueDict))

    def simpleBedUpdateWrapper(self, valueDict, db, genome, column):
        """
        If your details-mode classifier has BED records for values in its valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.updateRows(cur, genome, self.primaryKeyColumn, column, self.detailsEntryIter(valueDict.iteritems()))     

    def detailsEntryIter(self, valueIter):
        """
        General case for converting a valueIter to a string entry representing 1 or more BED records
        valueIters are generally lists of lists where each sublist represents a BED record
        """
        for aId, entry in valueIter:
            if entry is None:
                yield None, aId
            elif len(entry) == 0:
                raise RuntimeError("Empty list in details entry. This is not allowed.")
            elif type(entry[0]) != list:
                #only one entry
                yield "\t".join(map(str,entry)), aId
            else:
                bedEntries = ["\t".join(map(str, x)) for x in entry]
                yield "\n".join(bedEntries), aId

    def initializeDb(self, dbPath, classifiers, dataType=None):
        if dataType is None:
            columnDefinitions = [[x.__name__, x._getType()] for x in classifiers]
        else:
            columnDefinitions = [[x.__name__, dataType] for x in classifiers]
        #find alignment IDs from PSLs (primary key for database)
        for genome in self.genomes:
            aIds = set(x.split()[9] for x in open(self.alnPslDict[genome]))
            self.initializeSqlTable(dbPath, genome, columnDefinitions, self.primaryKeyColumn)
            self.initializeSqlRows(dbPath, genome, aIds, self.primaryKeyColumn)

    def initializeSqlTable(self, db, genome, columns, primaryKey):
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.initializeTable(cur, genome, columns, primaryKey)

    def initializeSqlRows(self, db, genome, aIds, primaryKey):
        with sql_lib.ExclusiveSqlConnection(db) as cur:
            sql_lib.insertRows(cur, genome, primaryKey, [primaryKey], izip_longest(aIds, [None]))