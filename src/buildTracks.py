import os
from itertools import izip_longest, product, izip
import sqlite3 as sql
from collections import defaultdict

import src.queries
import lib.sqlite_lib as sql_lib
import lib.psl_lib as psl_lib
from lib.general_lib import functionsInModule
from src.abstractClassifier import AbstractClassifier
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
    def __init__(self, outDir, genomes, primaryKeyColumn, sizes, beds, annotationBed):
        Target.__init__(self)
        self.outDir = outDir
        self.bedDir = os.path.join(self.outDir, "bedfiles")
        self.bigBedDir = os.path.join(self.outDir, "bigBedfiles")
        self.genomes = genomes
        self.primaryKeyColumn = primaryKeyColumn
        self.beds = beds
        self.sizes = sizes
        self.annotationBed = annotationBed
        self.categories = functionsInModule(src.queries)
        # bring in abstractClassifier colors
        self.colors = AbstractClassifier.colors

    def writeBed(self, genome, detailsFields, classifyFields, classifyValues, classifyOperations, categoryName):
        bedPath = os.path.join(self.bedDir, categoryName, genome, genome + ".bed")
        with open(bedPath, "w") as outf:
            for details in detailsFields:
                for record in sql_lib.selectBetweenDatabases(self.cur, "details", details, classifyFields, classifyValues, classifyOperations, self.primaryKeyColumn, genome):
                    if record[0] == None:
                        continue
                    elif type(record[0]) == type(u''):
                        outf.write(record[0] + "\n")
                    else:
                        for x in record[0]:
                            outf.write(x)+"\n"
        return bedPath

    def buildBigBed(self, bedPath, sizePath, genome, categoryName):
        bigBedPath = os.path.join(self.bigBedDir, categoryName, genome, genome + ".bb")
        system("bedSort {} {}".format(bedPath, bedPath))
        system("bedToBigBed -extraIndex=name {} {} {}".format(bedPath, sizePath, bigBedPath))

    def recolorTransMap(self, genome, bed):
        """
        Recolors the comparativeAnnotation results based on the scheme assembly > alignment > biology. Transcripts not in
        one of these categories will become black. Also removes the unique tag from transcript IDs.
        """
        records = [x.split() for x in open(bed)]
        # first we recolor everything black
        for x in records:
            x[8] = "0"
        # now we find all interesting biology and color that the interesting biology color
        detailsFields, classifyFields, classifyValues, classifyOperations = src.queries.interestingBiology()
        aIds = {x[0] for x in sql_lib.selectBetweenDatabases(self.cur, "details", self.primaryKeyColumn, classifyFields, classifyValues, classifyOperations, self.primaryKeyColumn, genome)}
        for x in records:
            if x[3] in aIds:
                x[8] = self.colors["mutation"]
        # now the alignment errors...
        detailsFields, classifyFields, classifyValues, classifyOperations = src.queries.alignmentErrors()
        aIds = {x[0] for x in sql_lib.selectBetweenDatabases(self.cur, "details", self.primaryKeyColumn, classifyFields, classifyValues, classifyOperations, self.primaryKeyColumn, genome)}
        for x in records:
            if x[3] in aIds:
                x[8] = self.colors["alignment"]
        # finally the assembly
        detailsFields, classifyFields, classifyValues, classifyOperations = src.queries.assemblyErrors()
        aIds = {x[0] for x in sql_lib.selectBetweenDatabases(self.cur, "details", self.primaryKeyColumn, classifyFields, classifyValues, classifyOperations, self.primaryKeyColumn, genome)}
        for x in records:
            if x[3] in aIds:
                x[8] = self.colors["assembly"]
        # remove alignment IDs
        for x in records:
            x[3] = psl_lib.removeAlignmentNumber(x[3])
        return ["\t".join(x) for x in records]

    def run(self):
        self.con = sql.connect(os.path.join(self.outDir, "classify.db"))
        self.cur = self.con.cursor()
        sql_lib.attachDatabase(self.con, os.path.join(self.outDir, "details.db"), "details")

        if not os.path.exists(self.bedDir):
            os.mkdir(self.bedDir)
        if not os.path.exists(self.bigBedDir):
            os.mkdir(self.bigBedDir)
        if not os.path.exists(os.path.join(self.bedDir, "comparativeAnnotation")):
            os.mkdir(os.path.join(self.bedDir, "comparativeAnnotation"))
        if not os.path.exists(os.path.join(self.bigBedDir, "comparativeAnnotation")):
            os.mkdir(os.path.join(self.bigBedDir, "comparativeAnnotation"))

        #build directory of comparativeAnnotation output
        for genome, bed, size in izip(self.genomes, self.beds, self.sizes):
            assert genome == os.path.basename(size).split(".")[0], (genome, os.path.basename(size).split(".")[0])
            if not os.path.exists(os.path.join(self.bedDir, "comparativeAnnotation", genome)):
                os.mkdir(os.path.join(self.bedDir, "comparativeAnnotation", genome))
            if not os.path.exists(os.path.join(self.bigBedDir, "comparativeAnnotation", genome)):
                os.mkdir(os.path.join(self.bigBedDir, "comparativeAnnotation", genome))
            recolored_records = self.recolorTransMap(genome, bed)
            new_bed_path = os.path.join(self.bedDir, "comparativeAnnotation", genome, genome + ".bed")
            with open(new_bed_path, 'w') as outf:
                for l in recolored_records:
                    outf.write(l + "\n")
            self.buildBigBed(new_bed_path, size, genome, "comparativeAnnotation")

        for category in self.categories:
            if not os.path.exists(os.path.join(self.bedDir, category.__name__)):
                os.mkdir(os.path.join(self.bedDir, category.__name__))
            if not os.path.exists(os.path.join(self.bigBedDir, category.__name__)):
                os.mkdir(os.path.join(self.bigBedDir, category.__name__))
            detailsFields, classifyFields, classifyValues, classifyOperations = category()
            for genome, sizePath in izip(self.genomes, self.sizes):
                if not os.path.exists(os.path.join(self.bedDir, category.__name__, genome)):
                    os.mkdir(os.path.join(self.bedDir, category.__name__, genome))
                if not os.path.exists(os.path.join(self.bigBedDir, category.__name__, genome)):
                    os.mkdir(os.path.join(self.bigBedDir, category.__name__, genome))
                bedPath = self.writeBed(genome, detailsFields, classifyFields, classifyValues, classifyOperations, category.__name__)
                # dumb - checks to make sure the BED is not empty so bedToBigBed doesn't crash
                if os.stat(bedPath).st_size != 0:
                    self.buildBigBed(bedPath, sizePath, genome, category.__name__)
