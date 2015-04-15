import os
import cPickle as pickle

from jobTree.src.bioio import logger
from src.abstractClassifier import Attribute

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import formatRatio

class TranscriptId(Attribute):
    """
    Creates a column representing the transcript Id
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        valueDict = {aId: psl_lib.removeAlignmentNumber(aId) for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneId(Attribute):
    """
    Creates a column representing the gene Id
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneId for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneName(Attribute):
    """
    Creates a column representing the gene name
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneName for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class GeneType(Attribute):
    """
    Creates a column representing the gene type
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneType for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class TranscriptType(Attribute):
    """
    Creates a column representing the transcript type
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAttributeDict()
        self.getAlignmentDict()
        valueDict = {aId: self.attributeDict[psl_lib.removeAlignmentNumber(aId)].transcriptType for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceChrom(Attribute):
    """
    Creates a column representing the source chromosome
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosome for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStart(Attribute):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].start for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStop(Attribute):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: self.annotationDict[psl_lib.removeAlignmentNumber(aId)].stop for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class SourceStrand(Attribute):
    """
    Creates a column representing the source genomic strand.
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        valueDict = {aId: seq_lib.convertStrand(self.annotationDict[psl_lib.removeAlignmentNumber(aId)].strand) for aId in self.alignmentDict}
        self.dumpValueDict(valueDict)


class DestChrom(Attribute):
    """
    Creates a column representing the dest chromosome
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].chromosome for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStart(Attribute):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].start for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStop(Attribute):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    @property
    def getType(self):
        return "INTEGER"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        valueDict = {aId: self.transcriptDict[aId].stop for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class DestStrand(Attribute):
    """
    Creates a column representing the dest genomic strand.
    """
    @property
    def getType(self):
        return "TEXT"

    def run(self):
        logger.info("Starting attribute {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        valueDict = {aId: seq_lib.convertStrand(self.transcriptDict[aId].strand) for aId in self.transcriptDict}
        self.dumpValueDict(valueDict)


class AlignmentCoverage(Attribute):
    """
    Calculates alignment coverage:

    (matches + mismatches) / qSize

    Reports the value as a REAL between 0 and 1
    """
    @property
    def getType(self):
        return "REAL"

    def run(self):
        logger.info("Starting attribute analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches + aln.misMatches, aln.qSize)
        self.dumpValueDict(valueDict)


class AlignmentIdentity(Attribute):
    """

    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1
    """
    @property
    def getType(self):
        return "REAL"

    def run(self):
        logger.info("Starting classifying analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches, aln.matches + aln.misMatches + aln.qNumInsert)
        self.dumpValueDict(valueDict)
