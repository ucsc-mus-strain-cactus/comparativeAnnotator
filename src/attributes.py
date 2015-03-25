import os
import cPickle as pickle

from jobTree.src.bioio import logger
from src.abstractClassifier import AbstractClassifier

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import formatRatio

class Attribute(AbstractClassifier):
    """Need to overwrite the dumpValueDict method for attributes"""
    def getAttributeDict(self):
        self.attributeDict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)        

    def dumpValueDict(self, valueDict):
        """
        Dumps a attribute dict.
        """
        #with open(os.path.join(self.getGlobalTempDir(), "Attribute" + self.genome), "wb") as outf:
        with open(os.path.join(self.outDir, "Attribute" + self.getColumn() + self.genome), "wb") as outf:
            pickle.dump(valueDict, outf)


class TranscriptId(Attribute):
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


class GeneId(Attribute):
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


class GeneName(Attribute):
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


class GeneType(Attribute):
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


class TranscriptType(Attribute):
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


class SourceChrom(Attribute):
    """
    Creates a column representing the source chromosome
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosome for aId in self.aIds}
        self.dumpValueDict(valueDict)


class SourceStart(Attribute):
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
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].start for aId in self.aIds}
        self.dumpValueDict(valueDict)


class SourceStop(Attribute):
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
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].stop for aId in self.aIds}
        self.dumpValueDict(valueDict)


class SourceStrand(Attribute):
    """
    Creates a column representing the source genomic strand.
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        valueDict = {aId: seq_lib.convertStrand(self.annotationDict[psl_lib.removeAlignmentNumber(aId)].strand)
                     for aId in self.aIds}
        self.dumpValueDict(valueDict)


class DestChrom(Attribute):
    """
    Creates a column representing the dest chromosome
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosome for aId in self.aIds if aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStart(Attribute):
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
        valueDict = {aId: self.transcriptDict[aId].start for aId in self.aIds if aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStop(Attribute):
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
        valueDict = {aId: self.transcriptDict[aId].stop for aId in self.aIds if aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStrand(Attribute):
    """
    Creates a column representing the dest genomic strand.
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {aId: seq_lib.convertStrand(self.transcriptDict[aId].strand) for aId in self.aIds if aId in 
                     self.transcriptDict}
        self.dumpValueDict(valueDict)

class AlignmentCoverage(Attribute):
    """
    Calculates alignment coverage:

    (matches + mismatches) / qSize

    Reports the value as a REAL between 0 and 1
    """
    @staticmethod
    def _getType():
        return "REAL"

    def run(self):
        logger.info("Starting attribute analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches + aln.misMatches, aln.qSize)
        logger.info(
            "Attribute {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class AlignmentIdentity(Attribute):
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
            "Attribute {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)

