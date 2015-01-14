from collections import defaultdict

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib


class AlignmentAbutsLeft(AbstractClassifier):
    """
    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """

    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):

        self.get_alignment_dict()

        s_dict = defaultdict(int)
        for a_id, aln in self.alignment_dict.iteritems():
            if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
                s_dict[a_id] = 1
            elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                s_dict[a_id] = 1

        self.upsert_dict_wrapper(s_dict)