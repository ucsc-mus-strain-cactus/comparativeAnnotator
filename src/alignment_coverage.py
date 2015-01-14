from src.abstract_classifier import AbstractClassifier
from lib.general_lib import formatRatio
import lib.sequence_lib as seq_lib

class AlignmentCoverage(AbstractClassifier):
    """

    Calculates alignment coverage:

    (matches + mismatches) / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    """

    @staticmethod
    def __type__():
        return "REAL"

    def run(self):
        self.get_alignment_dict()

        s_dict = {}
        for a_id, aln in self.alignment_dict.iteritems():
            r = formatRatio(aln.matches + aln.misMatches, aln.matches + aln.misMatches + aln.qNumInsert)
            assert r >= 0 and r <= 1
            s_dict[a_id] = r

        self.upsert_dict_wrapper(s_dict)