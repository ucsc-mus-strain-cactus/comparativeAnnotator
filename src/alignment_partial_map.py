from collections import defaultdict

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib


class AlignmentPartialMap(AbstractClassifier):
    """
    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """

    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_alignment_dict()

        s_dict = defaultdict(int)
        for a_id, aln in self.alignment_dict.iteritems():
            if aln.qSize != aln.qEnd - aln.qStart:
                s_dict[a_id] = 1

        self.upsert_dict_wrapper(s_dict)