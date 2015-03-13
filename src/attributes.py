from src.abstractClassifier import AbstractClassifier

from jobTree.src.bioio import logger

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import formatRatio


class TranscriptId(AbstractClassifier):
    """
    Creates a column representing the transcript Id
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        valueDict = {aId: psl_lib.removeAlignmentNumber(aId) for aId in self.aIds}
        self.dumpValueDict(valueDict)


class GeneId(AbstractClassifier):
    """
    Creates a column representing the gene Id
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAttributeDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneId for aId in self.aIds}
        self.dumpValueDict(valueDict)


class GeneName(AbstractClassifier):
    """
    Creates a column representing the gene name
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAttributeDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneName for aId in self.aIds}
        self.dumpValueDict(valueDict)


class GeneType(AbstractClassifier):
    """
    Creates a column representing the gene type
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAttributeDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneType for aId in self.aIds}
        self.dumpValueDict(valueDict)


class TranscriptType(AbstractClassifier):
    """
    Creates a column representing the transcript type
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAttributeDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].transcriptType for aId in self.aIds}
        self.dumpValueDict(valueDict)


class SourceChrom(AbstractClassifier):
    """
    Creates a column representing the source chromosome
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.chromosome for aId
                     in self.aIds}
        self.dumpValueDict(valueDict)


class SourceStart(AbstractClassifier):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.start for aId in
                     self.aIds}
        self.dumpValueDict(valueDict)


class SourceStop(AbstractClassifier):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.stop for aId in
                     self.aIds}
        self.dumpValueDict(valueDict)


class SourceStrand(AbstractClassifier):
    """
    Creates a column representing the source genomic strand.
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {
            aId: seq_lib.convertStrand(
                self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.strand)
            for aId in self.aIds}
        self.dumpValueDict(valueDict)


class DestChrom(AbstractClassifier):
    """
    Creates a column representing the dest chromosome
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosomeInterval.chromosome for aId in self.aIds if
                     aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStart(AbstractClassifier):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosomeInterval.start for aId in self.aIds if
                     aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStop(AbstractClassifier):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosomeInterval.stop for aId in self.aIds if
                     aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStrand(AbstractClassifier):
    """
    Creates a column representing the dest genomic strand.
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: seq_lib.convertStrand(self.transcriptDict[aId].chromosomeInterval.strand) for aId in self.aIds
                     if aId in self.transcriptDict}
        self.dumpValueDict(valueDict)

class AlignmentCoverage(AbstractClassifier):
    """

    Calculates alignment coverage:

    (matches + mismatches) / qSize

    Reports the value as a REAL between 0 and 1

    """

    @staticmethod
    def _getType():
        return "REAL"

    def run(self):
        logger.info("Starting classifying analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches + aln.misMatches, aln.qSize)
        logger.info(
            "Classify {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class AlignmentIdentity(AbstractClassifier):
    """

    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    """

    @staticmethod
    def _getType():
        return "REAL"

    def run(self):
        logger.info("Starting classifying analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches, aln.matches + aln.misMatches + aln.qNumInsert)
        logger.info(
            "Classify {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)

