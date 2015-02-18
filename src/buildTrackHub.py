import os
from itertools import izip_longest, product
import cPickle as pickle
import sqlite3 as sql

import src.classifiers, src.details, src.attributes
import lib.sqlite_lib as sql_lib
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger, system

class BuildTrackHub(Target):
    """
    Builds a track hub out of the databases. First, initializes a trackHub in the directory specified.
    Then, builds one huge bigBed file for each comparison designed.

    Each function that defines a comparison must create:
        a list of details fields to put in the BED record.
        a list of classify fields to filter on.
        a matching list of values that the classify field should have.
    """
    def __init__(self, outDir, trackHubDir, genomes, classifiers, details, attributes, primaryKeyColumn, trackHubName,
                 dataDir):
        Target.__init__(self)
        self.outDir = outDir
        self.trackHubDir = trackHubDir
        self.genomes = genomes
        self.classifiers = classifiers
        self.details = details
        self.attributes = attributes
        self.primaryKeyColumn = primaryKeyColumn
        self.trackHubName = trackHubName
        #self.tempDir = self.getGlobalTempDir()
        self.dataDir = dataDir
        #self.categories = [self.mutations, self.assemblyErrors]
        self.categories = [self.mutations]

    def mutations(self):
        detailsFields = ["CodingInsertions", "CodingDeletions", "CodingMult3Insertions", "CodingMult3Deletions",
                         "CdsNonCanonSplice", "UtrNonCanonSplice", "CdsUnknownSplice", "UtrUnknownSplice"]
        classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases",
                          "ScaffoldGap"]
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

    def startHub(self):
        if not os.path.exists(self.trackHubDir):
            os.mkdir(self.trackHubDir)
        with open(os.path.join(self.trackHubDir, "hub.txt"), "w") as outf:
            outf.write("hub {0}\nshortLabel {0}\nlongLabel {0}\ngenomesFile genomes.txt\nemail ian.t.fiddes@gmail.com\n"
                       .format(self.trackHubName))
        with open(os.path.join(self.trackHubDir, "genomes.txt"), "w") as outf:
            for genome in self.genomes:
                outf.write("genome {0}\ntrackDb {0}/trackDb.txt\ntwoBitPath\n{0}.2bit\n\n".format(genome))
            if not os.path.exists(os.path.join(self.trackHubDir, genome)):
                os.mkdir(os.path.join(self.trackHubDir, genome))
        for genome in self.genomes:
            if not os.path.exists(os.path.join(self.trackHubDir, genome)):
                os.mkdir(os.path.join(self.trackHubDir, genome))
            with open(os.path.join(self.trackHubDir, genome, genome + ".html"), "w") as outf:
                outf.write(genome + "\n")
        for genome in self.genomes:
            for ext in [".chrom.sizes", ".2bit"]:
                s = os.path.abspath(os.path.join(self.dataDir, genome + ext))
                t = os.path.join(self.trackHubDir, genome, genome + ext)
                if not os.path.exists(t):
                    system("ln -s {} {}".format(s, t))

    def writeBed(self, genome, detailsFields, classifyFields, classifyValues, classifyOperations, categoryName):
        bedPath = os.path.join(self.trackHubDir, genome, categoryName + ".bed")
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
            return bedPath

    def addTrack(self, name, genome, bigBedPath):
        with open(os.path.join(self.trackHubDir, genome, "trackDb.txt"), "a") as outf:
            outf.write("track {0}\nshortLabel {0} {1}\nlongLabel {0} {1}\nitemRgb on\ntype bigBed 12\nbigDataUrl {2}\n\n".format(genome,
                                                                                                                 name, os.path.basename(bigBedPath)))

    def buildBigBed(self, bedPath, name, genome):
        bigBedPath = os.path.join(self.trackHubDir, genome, name + ".bb")
        chromSizesPath = os.path.join(self.dataDir, genome + ".chrom.sizes")
        system("bedSort {0} {0}".format(bedPath))
        system("bedToBigBed {} {} {}".format(bedPath, chromSizesPath, bigBedPath))
        self.addTrack(name, genome, bigBedPath)

    def run(self):
        self.startHub()
        self.con = sql.connect(os.path.join(self.outDir, "classify.db"))
        self.cur = self.con.cursor()
        sql_lib.attachDatabase(self.con, os.path.join(self.outDir, "details.db"), "details")
        for category in self.categories:
            detailsFields, classifyFields, classifyValues, classifyOperations = category()
            for genome in self.genomes:
                bedPath = self.writeBed(genome, detailsFields, classifyFields, classifyValues, classifyOperations,
                                        category.__name__)
                self.buildBigBed(bedPath, category.__name__, genome)


