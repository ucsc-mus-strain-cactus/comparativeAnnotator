from collections import Counter
from itertools import izip

from src.abstract_classifier import AbstractClassifier
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class UnknownBases(AbstractClassifier):
    """
    Counts the number of Ns in the target sequence within alignment blocks

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        #call methods from abstract_classifier that pull data needed for this classifier
        self.get_alignments()
        self.get_seq_dict()

        counts = Counter()
        for aln in self.alignments:
            if aln.strand == "+":
                for tStart, blockSize in izip(aln.tStarts, aln.blockSizes):
                    seq = self.seq_dict[aln.tName][tStart : tStart + blockSize]
                    counts[aln.qName] += seq.count("N")
            else:
                #on negative strand the tStarts are (+) strand but blockSizes are in
                #transcript orientation
                for tStart, blockSize in izip(aln.tStarts, reversed(aln.blockSizes)):
                    seq = self.seq_dict[aln.tName][tStart : tStart + blockSize]
                    counts[aln.qName] += seq.count("N")

        self.upsert_dict_wrapper(counts)