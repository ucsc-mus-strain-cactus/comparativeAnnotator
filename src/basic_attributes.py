from src.abstract_classifier import AbstractClassifier
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class TranscriptID(AbstractClassifier):
    """
    Creates a column representing the transcript ID
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):
        #call methods from abstract_classifier that pull data needed for this classifier
        self.get_alignments()

        t_map = {}
        for aln in self.alignments:
            t_map[aln.qName] = psl_lib.removeAlignmentNumber(aln.qName)

        self.upsert_dict_wrapper(t_map)


class GeneID(AbstractClassifier):
    """
    Creates a column representing the gene ID
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):
        self.get_alignments()
        self.get_transcript_attributes()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            t_map[aln.qName] = self.attribute_dict[tName].geneID

        self.upsert_dict_wrapper(t_map)


class GeneName(AbstractClassifier):
    """
    Creates a column representing the gene name
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):
        self.get_alignments()
        self.get_transcript_attributes()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            t_map[aln.qName] = self.attribute_dict[tName].geneName

        self.upsert_dict_wrapper(t_map)


class GeneType(AbstractClassifier):
    """
    Creates a column representing the gene type
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):
        self.get_alignments()
        self.get_transcript_attributes()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            t_map[aln.qName] = self.attribute_dict[tName].geneType

        self.upsert_dict_wrapper(t_map)    


class TranscriptType(AbstractClassifier):
    """
    Creates a column representing the transcript type
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):
        self.get_alignments()
        self.get_transcript_attributes()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            t_map[aln.qName] = self.attribute_dict[tName].transcriptType

        self.upsert_dict_wrapper(t_map)