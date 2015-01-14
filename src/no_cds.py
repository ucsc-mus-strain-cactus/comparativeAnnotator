from collections import defaultdict

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class NoCds(AbstractClassifier):
    """

    Looks to see if this transcript actually has a CDS

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_transcript_dict()

        s_dict = defaultdict(int)
        for a, t in self.transcript_dict.iteritems():
            if t.getCdsLength() < 3:
                s_dict[a] = 1

        self.upsert_dict_wrapper(s_dict)