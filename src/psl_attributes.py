from src.abstract_classifier import AbstractClassifier
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class SourceChrom(AbstractClassifier):
    """
    Creates a column representing the source chromosome
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_alignments()
        self.get_original_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            if tName not in self.original_transcript_dict:
                chrom = None
            else:
                chrom = self.original_transcript_dict[tName].chromosomeInterval.chromosome
            t_map[aln.qName] = chrom

        self.upsert_dict_wrapper(t_map)


class SourceStart(AbstractClassifier):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):   
        self.get_alignments()
        self.get_original_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            if tName not in self.original_transcript_dict:
                pos = None
            else:
                pos = self.original_transcript_dict[tName].chromosomeInterval.start
            t_map[aln.qName] = pos

        self.upsert_dict_wrapper(t_map)    


class SourceStop(AbstractClassifier):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):   
        self.get_alignments()
        self.get_original_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            if tName not in self.original_transcript_dict:
                pos = None
            else:
                pos = self.original_transcript_dict[tName].chromosomeInterval.stop
            t_map[aln.qName] = pos

        self.upsert_dict_wrapper(t_map)    


class SourceStrand(AbstractClassifier):
    """
    Creates a column representing the source genomic strand.
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):   
        self.get_alignments()
        self.get_original_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            tName = psl_lib.removeAlignmentNumber(aln.qName)
            if tName not in self.original_transcript_dict:
                s = None
            else:
                s = self.original_transcript_dict[tName].chromosomeInterval.strand
            if s is True:
                t_map[aln.qName] = "+"
            elif s is False:
                t_map[aln.qName] = "-"

        self.upsert_dict_wrapper(t_map)    


class DestChrom(AbstractClassifier):
    """
    Creates a column representing the dest chromosome
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):
        self.get_alignments()
        self.get_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            if aln.qName not in self.transcript_dict:
                chrom = None
            else:
                chrom = self.transcript_dict[aln.qName].chromosomeInterval.chromosome
            t_map[aln.qName] = chrom

        self.upsert_dict_wrapper(t_map)


class DestStart(AbstractClassifier):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):   
        self.get_alignments()
        self.get_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            if aln.qName not in self.transcript_dict:
                pos = None
            else:
                pos = self.transcript_dict[aln.qName].chromosomeInterval.start
            t_map[aln.qName] = pos

        self.upsert_dict_wrapper(t_map)    


class DestStop(AbstractClassifier):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def run(self):   
        self.get_alignments()
        self.get_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            if aln.qName not in self.transcript_dict:
                pos = None
            else:
                pos = self.transcript_dict[aln.qName].chromosomeInterval.stop
            t_map[aln.qName] = pos

        self.upsert_dict_wrapper(t_map)    


class DestStrand(AbstractClassifier):
    """
    Creates a column representing the dest genomic strand.
    """
    @staticmethod
    def __type__():
        return "TEXT"

    def run(self):   
        self.get_alignments()
        self.get_transcript_dict()

        t_map = {}
        for aln in self.alignments:
            if aln.qName not in self.transcript_dict:
                s = None
            else:
                s = self.transcript_dict[aln.qName].chromosomeInterval.strand
            if s is True:
                t_map[aln.qName] = "+"
            elif s is False:
                t_map[aln.qName] = "-"

        self.upsert_dict_wrapper(t_map)            