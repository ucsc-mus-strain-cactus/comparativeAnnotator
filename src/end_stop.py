from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class EndStop(AbstractClassifier):
    """
    Looks at the end of the coding region (thickEnd) and sees if the last
    three bases are a stop codon ('TAA', 'TGA', 'TAG')

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
            s = t.getProteinSequence(self.seq_dict)
            if len(s) > 0 and s[-1] != "*":
                s_dict[a] = 1
            else:
                s_dict[a] = 0

        self.upsert_dict_wrapper(s_dict)