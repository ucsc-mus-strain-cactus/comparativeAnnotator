from collections import defaultdict

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class CdsMult3Gap(AbstractClassifier):
    """

    Are any of the short CDS introns a multiple of 3?

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self, short_intron_size=30):
        self.get_transcript_dict()

        s_dict = defaultdict(int)
        for a, t in self.transcript_dict.iteritems():
            for i in xrange(len(t.intronIntervals)):
                if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                    if len(t.intronIntervals[i]) <= short_intron_size and len(t.intronIntervals[i]) % 3 == 0:
                        s_dict[a] = 1
                        break

        self.upsert_dict_wrapper(s_dict)