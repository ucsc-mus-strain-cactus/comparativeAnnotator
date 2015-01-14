from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib

class BeginStart(AbstractClassifier):
    """

    Reports the 0-based position of the start codon in the transcript.

    Writes -1 to the database if this transcript either:
    1) has no thick region (thickStart == thickStop)
    2) has a thick region smaller than 3
    3) has no start codon in the thick window

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_transcript_dict()
        self.get_seq_dict()

        s_dict = {}
        for a, t in self.transcript_dict.iteritems():
            s = t.getCds(self.seq_dict)
            #ATG is the only start codon
            if len(s) == 0 or s[:3] != "ATG":
                s_dict[a] = -1
            else:
                s_dict[a] = t.cdsCoordinateToTranscript(0)

        self.upsert_dict_wrapper(s_dict)