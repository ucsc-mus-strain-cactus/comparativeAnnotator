import os
from itertools import izip_longest, product
import sqlite3 as sql
from collections import defaultdict

from queries import *
from src.summaryStatistics import SummaryStatistics

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
        a matching list of AND/OR operations to link the logic together
    """
    def __init__(self, outDir, genomes, primaryKeyColumn, dataDir, geneCheckBedDict, annotationBed):
        Target.__init__(self)
        self.outDir = outDir
        self.bedDir = os.path.join(self.outDir, "bedfiles")
        self.genomes = genomes
        self.primaryKeyColumn = primaryKeyColumn
        self.geneCheckBedDict = geneCheckBedDict
        self.dataDir = dataDir
        self.annotationBed = annotationBed
        self.categories = [mutations, inFrameStop, interestingBiology, assemblyErrors, alignmentErrors, everything]

    def writeBed(self, genome, detailsFields, classifyFields, classifyValues, classifyOperations, categoryName):
        bedPath = os.path.join(self.bedDir, categoryName, genome, genome + ".bed")
        with open(bedPath, "w") as outf:
            for details in detailsFields:
                for record in sql_lib.selectBetweenDatabases(self.cur, "details", details, classifyFields, 
                                                             classifyValues, classifyOperations, self.primaryKeyColumn, 
                                                             genome):
                    if record[0] == None:
                        continue
                    elif type(record[0]) == type(u''):
                        outf.write(record[0] + "\n")
                    else:
                        for x in record[0]:
                            outf.write(x)+"\n"
        return bedPath

    def buildBigBed(self, bedPath, genome, categoryName):
        bigBedPath = os.path.join(self.bedDir, categoryName, genome, genome + ".bb")
        chromSizesPath = os.path.join(self.dataDir, genome + ".chrom.sizes")
        #system("bedSort {} {}".format(bedPath, os.path.join(self.getLocalTempDir(), "tmp"))
        system("bedSort {} {}".format(bedPath, "tmp"))
        #system("bedToBigBed {} {} {}".format(os.path.join(self.getLocalTempDir(), "tmp"), chromSizesPath, bigBedPath))
        system("bedToBigBed {} {} {}".format("tmp", chromSizesPath, bigBedPath))
        os.remove("tmp")

    def run(self):
        if not os.path.exists(self.bedDir):
            os.mkdir(self.bedDir)
        if not os.path.exists(os.path.join(self.bedDir, "transMap")):
            os.mkdir(os.path.join(self.bedDir, "transMap"))
        
        #build directory of transMap output
        for genome, bed in self.geneCheckBedDict.iteritems():
            if not os.path.exists(os.path.join(self.bedDir, "transMap", genome)):
                os.mkdir(os.path.join(self.bedDir, "transMap", genome))
            self.buildBigBed(bed, genome, "transMap")
        #don't need if including self in original analysis
        if not os.path.exists(os.path.join(self.bedDir, "transMap", "C57B6J")):
            os.mkdir(os.path.join(self.bedDir, "transMap", "C57B6J"))
        if not os.path.exists(os.path.join(self.bedDir, "transMap", "C57B6J", "C57B6J" + ".bed")):
            system("ln -s {} {}".format(os.path.abspath(self.annotationBed), os.path.join(self.bedDir, "transMap", "C57B6J", "C57B6J" + ".bed")))

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
                bedPath = self.writeBed(genome, detailsFields, classifyFields, classifyValues, classifyOperations,
                                        category.__name__)
                #dumb
                if len(open(bedPath).readlines()) > 0:
                    self.buildBigBed(bedPath, genome, category.__name__)
                os.remove(bedPath)

        self.setFollowOnTarget(SummaryStatistics(self.outDir, self.genomes))


