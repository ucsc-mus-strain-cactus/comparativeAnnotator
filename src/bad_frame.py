from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class BadFrame(AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_transcript_dict()
        self.get_seq_dict()

        s_dict = {}
        for a, t in self.transcript_dict.iteritems():
            if t.getCdsLength() % 3 != 0:
                s_dict[a] = 1
            else:
                s_dict[a] = 0

        self.upsert_dict_wrapper(s_dict)