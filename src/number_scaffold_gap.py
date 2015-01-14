from collections import defaultdict
import re

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class NumberScaffoldGap(AbstractClassifier):
    """

    How many 100bp N runs are there in the alignment?
    100bp N runs are markers of scaffold gaps.

    """

    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_alignment_dict()
        self.get_seq_dict()

        r = re.compile("[N]{100}")

        s_dict = defaultdict(int)
        for a_id, aln in self.alignment_dict.iteritems():
            dest_seq = self.seq_dict[aln.tName][aln.tStart : aln.tEnd].upper()
            if re.search(r, dest_seq) is not None:
                s_dict[a_id] = 1

        self.upsert_dict_wrapper(s_dict)