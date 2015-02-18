import os
from itertools import izip_longest, product
import cPickle as pickle
import sqlite3 as sql

import src.classifiers, src.details, src.attributes
import lib.sqlite_lib as sql_lib
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger, system

class BuildTracks(Target):
    """
    Builds a track hub out of the databases. First, initializes a trackHub in the directory specified.
    Then, builds one huge bigBed file for each comparison designed.

    Each function that defines a comparison must create:
        a list of details fields to put in the BED record.
        a list of classify fields to filter on.
        a matching list of values that the classify field should have.
    """
    def __init__(self, outDir, trackHubDir, genomes, classifiers, details, attributes, primaryKeyColumn,
                      dataDir, geneCheckBedDict):
        Target.__init__(self)
        self.outDir = outDir
        self.trackHubDir = trackHubDir
        self.bedDir = os.path.join(self.trackHubDir, "bedfiles")
        self.genomes = genomes
        self.classifiers = classifiers
        self.details = details
        self.attributes = attributes
        self.primaryKeyColumn = primaryKeyColumn
        self.geneCheckBedDict = geneCheckBedDict
        #self.tempDir = self.getGlobalTempDir()
        self.dataDir = dataDir
        #self.categories = [self.mutations, self.assemblyErrors]
        self.categories = [self.mutations]

    def mutations(self):
        detailsFields = ["CodingInsertions", "CodingDeletions", "CodingMult3Insertions", "CodingMult3Deletions",
                         "CdsNonCanonSplice", "UtrNonCanonSplice", "CdsUnknownSplice", "UtrUnknownSplice", "CdsMult3Gap"]
        classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases",
                          "ScaffoldGap", "UtrGap", "CdsGap"]
        classifyOperations = ["AND"] * len(classifyFields)
        classifyValues = [0] * len(classifyFields)
        return detailsFields, classifyFields, classifyValues, classifyOperations

    def assemblyErrors(self):
        detailsFields = [x.__name__ for x in self.details]
        classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases",
                          "ScaffoldGap"]
        classifyOperations = ["OR"] * len(classifyFields)
        classifyValues = [1] * len(classifyFields)
        return detailsFields, classifyFields, classifyValues, classifyOperations

    def writeBed(self, genome, detailsFields, classifyFields, classifyValues, classifyOperations, categoryName):
        bedPath = os.path.join(self.bedDir, categoryName, genome, genome + ".bed")
        with open(bedPath, "w") as outf:
            for details in detailsFields:
                for record in sql_lib.selectBetweenDatabases(self.cur, "details", details, classifyFields, classifyValues,
                                                              classifyOperations, self.primaryKeyColumn, genome):
                    if record[0] == None:
                        continue
                    elif type(record[0]) == type(u''):
                        outf.write(record[0] + "\n")
                    else:
                        for x in record[0]:
                            outf.write(x)+"\n"
        self.buildBigBed(bedPath, os.path.join(self.bedDir, categoryName, genome, genome + ".bb"), genome)
        os.remove(bedPath)

    def buildBigBed(self, bedPath, bigBedPath, genome):
        chromSizesPath = os.path.join(self.dataDir, genome + ".chrom.sizes")
        #tmp = os.path.join(self.tempDir, "tmp")
        tmp = os.path.join(self.bedDir, "tmp")
        system("bedSort {} {}".format(bedPath, tmp))
        system("bedToBigBed {} {} {}".format(tmp, chromSizesPath, bigBedPath))
        os.remove(tmp)

    def run(self):
        if not os.path.exists(self.bedDir):
            os.mkdir(self.bedDir)
        if not os.path.exists(os.path.join(self.bedDir, "geneCheck")):
            os.mkdir(os.path.join(self.bedDir, "geneCheck"))
        #build directory of geneCheck output
        for genome, bed in self.geneCheckBedDict.iteritems():
            self.buildBigBed(bed, os.path.join(self.bedDir, "geneCheck", genome + ".bb"), genome)
        self.con = sql.connect(os.path.join(self.outDir, "classify.db"))
        self.cur = self.con.cursor()
        sql_lib.attachDatabase(self.con, os.path.join(self.outDir, "details.db"), "details")
        for category in self.categories:
            if not os.path.exists(os.path.join(self.bedDir, category.__name__)):
                os.mkdir(os.path.join(self.bedDir, category.__name__))
            detailsFields, classifyFields, classifyValues, classifyOperations = category()
            for genome in self.genomes:
                if not os.path.exists(os.path.join(self.bedDir, category.__name__, genome)):
                    os.mkdir(os.path.join(self.bedDir, category.__name__, genome))
                self.writeBed(genome, detailsFields, classifyFields, classifyValues, classifyOperations,
                                        category.__name__)


