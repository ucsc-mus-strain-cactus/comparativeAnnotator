import os
import sqlite3 as sql
from collections import defaultdict, Counter
from itertools import izip

import src.queries
from lib.general_lib import functionsInModule

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

        self.statistics = defaultdict(list)

        categories = functionsInModule(src.queries)

        for genome in self.genomes:
            self.statistics["coverage"].append(self.coverage(self.cur, genome))
            self.statistics["identity"].append(self.identity(self.cur, genome))
            self.statistics["paralogy"].append(self.paralogy(self.cur, genome))
            for category in categories:
                detailsFields, classifyFields, classifyValues, classifyOperations = category()
                self.statistics[category.__name__].append(round(100.0 * self.numberCategorized(self.cur, genome, classifyFields, detailsFields, classifyValues, classifyOperations) / self.numberRows(self.cur, genome), 3))

        with open(os.path.join(self.outDir, self.outDir.replace('/','') + "summary.tsv"), "w") as outf:
            outf.write("genomes\t"+"\t".join(self.genomes)+"\n")
            for category in ["coverage", "identity", "paralogy"] + [x.__name__ for x in categories]:
                outf.write(category + "\t" + "\t".join(map(str, self.statistics[category])) + "\n")


    def paralogy(self, cur, genome):
        cmd = """SELECT attributes.'{0}'.GeneName FROM attributes.'{0}'""".format(genome)
        genes = Counter([x[0] for x in cur.execute(cmd).fetchall()])
        return round(100.0 * len([x for x in genes if genes[x] > 1]) / len(genes), 3)

    def coverage(self, cur, genome):
        """
        Finds the best coverage percent for each gene and records the average.
        TODO: histogram
        """
        cmd = """SELECT main.'{0}'.AlignmentCoverage, attributes.'{0}'.GeneName FROM attributes.'{0}' JOIN main.'{0}' USING ('AlignmentId') WHERE main.'{0}'.AlignmentCoverage IS NOT NULL""".format(genome)
        r = cur.execute(cmd).fetchall()
        t = defaultdict(float)
        for val, gene in r:
            if t[gene] < val:
                t[gene] = val
        return round(100.0 * sum(t.values())/len(t), 3)

    def identity(self, cur, genome):
        """
        Finds the best identity percent for each gene and records the average.
        TODO: histogram
        """
        cmd = """SELECT main.'{0}'.AlignmentIdentity, attributes.'{0}'.GeneName FROM attributes.'{0}' JOIN main.'{0}' USING ('AlignmentId') WHERE main.'{0}'.AlignmentIdentity IS NOT NULL""".format(genome)
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