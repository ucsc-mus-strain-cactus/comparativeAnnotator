from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class InFrameStop(AbstractClassifier):
    """
    Looks for stop codons in frame in the coding sequence.

    Reports on in frame stop codons for each transcript.

    Records the 0-based transcript coordinate of an in-frame stop codon
    if it exists. Otherwise, records -1.

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_transcript_dict()
        self.get_seq_dict()

        s_dict = {}
        for a, t in self.transcript_dict.iteritems():
            #make sure this transcript has CDS
            #and more than 2 codons - can't have in frame stop without that
            cds_size = t.getCdsLength()
            if cds_size >= 9:
                for i in xrange(3, cds_size - 3, 3):
                    c = t.cdsCoordinateToAminoAcid(i, self.seq_dict)
                    if c == "*":
                        s_dict[a] = i
            else:
                s_dict[a] = -1

        self.upsert_dict_wrapper(s_dict)