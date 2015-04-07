from collections import defaultdict

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

from jobTree.src.bioio import logger
from src.abstractClassifier import AbstractClassifier

class GapFinder(AbstractClassifier):
    """
    Finds improperly sized introns.
    """

    def gapFinder(self, t, shortIntronSize, mult3, coding):
        """
        Looks for introns that are smaller than shortIntronSize.
        Reports on these depending on the mult3 and coding flags
        Introns that contain unknown bases will not be counted as they are likely alignment issues.
        TODO: make this its own classifier/flag for this classifier and call it 'UnknownGap'
        mult3: True for only mult3, False for no mult3, None for either.
        coding: True for only coding, False for non-coding, None for either
        """
        records = []
        for i, intron in enumerate(t.intronIntervals):
            if len(intron) >= shortIntronSize:
                continue
            if "N" in intron.getSequence(self.seqDict):
                continue
            if (intron.start >= t.thickStart and intron.stop <= t.thickStop) and (coding is True or coding is None):
                if len(intron) % 3 == 0 and (mult3 is True or mult3 is None):
                    records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.getColumn()))
                elif len(intron) % 3 != 0 and (mult3 is False or mult3 is None):
                    records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.getColumn()))
            elif ((intron.stop <= t.thickStart and intron.start <= t.thickStart) or (intron.start >= t.thickStop and \
                        intron.stop >= t.thickStop)) and (coding is False or coding is None):
                if len(intron) % 3 == 0 and (mult3 is True or mult3 is None):
                    records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.getColumn()))
                elif len(intron) % 3 != 0 and (mult3 is False or mult3 is None):
                    records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.getColumn()))
        return records

    def run(self, shortIntronSize=30, mult3=None, coding=None):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            records = self.gapFinder(t, shortIntronSize, mult3, coding)
            if len(records) == 0:
                classifyDict[aId] = 0
            else:
                classifyDict[aId] = 1
                detailsDict[aId] = records
        self.dumpValueDicts(classifyDict, detailsDict)