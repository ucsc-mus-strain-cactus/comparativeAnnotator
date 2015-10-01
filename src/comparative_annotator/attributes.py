from src.abstractClassifier import Attribute

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import format_ratio


class TranscriptId(Attribute):
    """
    Creates a column representing the transcript Id
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAlignmentDict()
        valueDict = {aId: psl_lib.remove_alignment_number(aId) for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneId(Attribute):
    """
    Creates a column representing the gene Id
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.remove_alignment_number(aId)].geneId for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneName(Attribute):
    """
    Creates a column representing the gene name
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.remove_alignment_number(aId)].geneName for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneType(Attribute):
    """
    Creates a column representing the gene type
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.remove_alignment_number(aId)].geneType for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class TranscriptType(Attribute):
    """
    Creates a column representing the transcript type
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.remove_alignment_number(aId)].transcriptType for aId in
                     self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceChrom(Attribute):
    """
    Creates a column representing the source chromosome
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.remove_alignment_number(aId)].chromosome for aId in
                     self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStart(Attribute):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.remove_alignment_number(aId)].start for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStop(Attribute):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.remove_alignment_number(aId)].stop for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStrand(Attribute):
    """
    Creates a column representing the source genomic strand.
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: seq_lib.convert_strand(self.annotationDict[psl_lib.remove_alignment_number(aId)].strand) for
                     aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class DestChrom(Attribute):
    """
    Creates a column representing the dest chromosome
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosome for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStart(Attribute):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].start for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStop(Attribute):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    @staticmethod
    def dataType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].stop for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStrand(Attribute):
    """
    Creates a column representing the dest genomic strand.
    """
    @staticmethod
    def dataType():
        return "TEXT"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId: seq_lib.convert_strand(self.transcriptDict[aId].strand) for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class AlignmentCoverage(Attribute):
    """
    Calculates alignment coverage:

    (matches + mismatches + repeat matches) / qSize

    Reports the value as a REAL between 0 and 1
    """
    @staticmethod
    def dataType():
        return "REAL"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = format_ratio(aln.matches + aln.misMatches + aln.repMatches, aln.qSize)
        self.dumpValueDict(valueDict)


class AlignmentIdentity(Attribute):
    """
    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1
    """
    @staticmethod
    def dataType():
        return "REAL"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = format_ratio(aln.matches + aln.repMatches, aln.matches + aln.repMatches + aln.misMatches +
                                         aln.qNumInsert)
        self.dumpValueDict(valueDict)
