import os
import sqlite3 as sql
from collections import defaultdict, Counter
from itertools import izip, izip_longest

import src.queries
from lib.general_lib import functionsInModule, DefaultOrderedDict

import lib.sqlite_lib as sql_lib

from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger, system

class SummaryStatistics(Target):
    def __init__(self, outDir, genomes):
        Target.__init__(self)
        self.outDir = outDir
        self.genomes = sorted(genomes)

    def run(self):
        self.con = sql.connect(os.path.join(self.outDir, "classify.db"))
        self.cur = self.con.cursor()

        sql_lib.attachDatabase(self.con, os.path.join(self.outDir, "details.db"), "details")
        sql_lib.attachDatabase(self.con, os.path.join(self.outDir, "attributes.db"), "attributes")

        self.statistics = DefaultOrderedDict(list)

        categories = functionsInModule(src.queries)

        for genome in self.genomes:
            self.statistics["identity"].append(self.identity(self.cur, genome))
            self.statistics["paralogy"].append(self.paralogy(self.cur, genome))
            self.statistics["coding_OK"].append(self.ok_coding(self.cur, genome))
            for biotype in ["protein_coding", "miRNA", "lincRNA", "processed_transcript"]:
                self.statistics[biotype + "_coverage"].append(self.coverage_biotype(self.cur, genome, biotype))
            for category in categories:
                detailsFields, classifyFields, classifyValues, classifyOperations = category()
                self.statistics[category.__name__].append(round(100.0 * self.numberCategorized(self.cur, genome, classifyFields, detailsFields, classifyValues, classifyOperations) / self.numberRows(self.cur, genome), 3))

        with open(os.path.join(self.outDir, self.outDir.split("/")[-1] + "summary.tsv"), "w") as outf:
            outf.write("genomes\t"+"\t".join(self.genomes)+"\n")
            for category in self.statistics:
                outf.write(category + "\t" + "\t".join(map(str, self.statistics[category])) + "\n")

    def paralogy(self, cur, genome):
        """
        Finds the number of paralogous alignments. This is defined as the number of gene IDs with more than one
        alignment.
        """
        cmd = """SELECT attributes.'{0}'.TranscriptId FROM attributes.'{0}'""".format(genome)
        genes = Counter([x[0] for x in cur.execute(cmd).fetchall()])
        return round(100.0 * len([x for x in genes if genes[x] > 1]) / len(genes), 3)

    def coverage_biotype(self, cur, genome, biotype):
        """
        Finds the best coverage percent for each TRANSCRIPT and records the average.
        """
        cmd = """SELECT attributes.'{0}'.AlignmentCoverage, attributes.'{0}'.TranscriptId FROM attributes.'{0}' WHERE attributes.'{0}'.'GeneType' = '{1}'""".format(genome, biotype)
        r = cur.execute(cmd).fetchall()
        t = defaultdict(float)
        for val, gene in r:
            if t[gene] < val:
                t[gene] = val
        return round(100.0 * sum(t.values())/len(t), 3)

    def coverage_biotypes(self, cur, genome, biotypes):
        cmd = """SELECT attributes.'{0}'.AlignmentCoverage, attributes.'{0}'.TranscriptId FROM attributes.'{0}' WHERE attributes.'{0}'.'GeneType' = '{1}'""".format(genome, biotype)
        t = defaultdict(float)
        for biotype in biotypes:
            r = cur.execute(cmd).fetchall()
            for val, gene in r:
                if t[gene] < val:
                    t[gene] = val
        return round(100.0 * sum(t.values())/len(t), 3)

    def identity(self, cur, genome):
        """
        Finds the best identity percent for each gene and records the average.
        TODO: histogram
        """
        cmd = """SELECT attributes.'{0}'.AlignmentIdentity, attributes.'{0}'.GeneName FROM attributes.'{0}' WHERE attributes.'{0}'.AlignmentIdentity IS NOT NULL""".format(genome)
        r = cur.execute(cmd).fetchall()
        t = defaultdict(float)
        for val, gene in r:
            if t[gene] < val:
                t[gene] = val
        return round(100.0 * sum(t.values())/len(t), 3)

    def numberCategorized(self, cur, genome, classifyFields, detailsFields, classifyValues, classifyOperations):
        """
        Finds the number of ALIGNMENTS categorized by a categorizing function
        """
        cmd = """SELECT COUNT(*) FROM main.'{0}' WHERE (""".format(genome)
        for col, mod in izip(classifyFields[:-1], classifyOperations):
            cmd += " main.'{}'.'{}' = ? {}".format(genome, col, mod)
        cmd += " main.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
        if len(detailsFields) == 1:
            cmd += " AND main.'{}'.'{}' = 1".format(genome, detailsFields[0])
        q = cur.execute(cmd, classifyValues)
        return int(q.fetchone()[0])

    def numberRows(self, cur, genome):
        return int(cur.execute("SELECT COUNT(*) FROM main.'{}'".format(genome)).fetchone()[0])

    def ok_coding(self, cur, genome):
        """
        Finds the number of TRANSCRIPTS whose best alignment has no problems.
        
        OK is defined as transcripts who do not hit any classifier except for nonsynonymous/synonymous/noncanon splices
        """
        biotype = "protein_coding"
        classifyFields = ["CodingInsertions","CodingDeletions", "CodingDeletions", "StartOutOfFrame", "FrameShift", "AlignmentAbutsLeft", "AlignmentAbutsRight", 
                          "AlignmentPartialMap", "BadFrame", "BeginStart", "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                          "EndStop", "InFrameStop", "ShortCds", "UnknownBases"]
        cmd = """SELECT main.'{0}'.'AlignmentId' FROM main.'{0}' WHERE (""".format(genome)
        for col in classifyFields[:-1]:
            cmd += " main.'{}'.'{}' = ? {}".format(genome, col, "AND")
        cmd += " main.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
        vals = [0] * len(classifyFields)
        # find the best alignment for each transcript
        b_cmd = """SELECT attributes.'{0}'.AlignmentCoverage, attributes.'{0}'.TranscriptId, attributes.'{0}'.AlignmentId FROM attributes.'{0}' WHERE attributes.'{0}'.GeneType = '{1}'""".format(genome, biotype)
        r = cur.execute(b_cmd).fetchall()
        ids = defaultdict(str)
        cov = defaultdict(float)
        for val, gene, alignment_id in r:
            if cov[gene] < val:
                cov[gene] = val
                ids[gene] = alignment_id
        ids = {y: x for x, y in ids.iteritems()}
        ok_alignments = cur.execute(cmd, vals).fetchall()
        ok_transcripts = [ids[x[0]] for x in ok_alignments if x[0] in ids]
        return round(100.0 * len(ok_transcripts) / len(ids), 3)