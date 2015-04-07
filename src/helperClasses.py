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


class SpliceSiteAnalysis(AbstractClassifier):
    """
    Analyzes intronIntervals for abnormal splice sites.
    If canonical is True, reports on all splice sites that are not GT:AG
    If canonical is False, reports on splice sites that do not fit any of the known splice sites.

    coding: True for only coding, False for non-coding, None for either

    This classifier is only applied to introns which are longer than a minimum intron size.
    """
    canonical = {"GT": "AG"}
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}

    def canonicalSplice(self, donor, acceptor):
        """Is this splice site canonical?"""
        if donor in self.canonical and self.canonical[donor] != acceptor:
            return False
        elif donor not in self.canonical:
            return False
        else:
            return True

    def unknownSplice(self, donor, acceptor):
        """Is this splice site neither canonical or noncanonical?"""
        if donor in self.non_canonical and self.non_canonical[donor] != acceptor:
            return False
        elif donor not in self.non_canonical:
            return False
        else:
            return True

    def run(self, canonical=False, coding=None, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for intron in t.intronIntervals:
                if len(intron) <= shortIntronSize:
                    continue
                if (intron.start >= t.thickStart and intron.stop <= t.thickStop) and (coding is True or coding is None):
                    seq = intron.getSequence(self.seqDict, strand=True)
                    if canonical is True and self.canonicalSplice(seq[:2], seq[-2:]) is False:
                        classifyDict[aId] = 1
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb(), self.getColumn()))
                    elif canonical is False and self.unknownSplice(seq[:2], seq[-2:]) is False:
                        classifyDict[aId] = 1
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb(), self.getColumn()))
                elif ((intron.stop <= t.thickStart and intron.start <= t.thickStart) or (intron.start >= t.thickStop and \
                        intron.stop >= t.thickStop)) and (coding is False or coding is None):
                    seq = intron.getSequence(self.seqDict, strand=True)
                    if canonical is True and self.canonicalSplice(seq[:2], seq[-2:]) is False:
                        classifyDict[aId] = 1
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb(), self.getColumn()))
                    elif canonical is False and self.unknownSplice(seq[:2], seq[-2:]) is False:
                        classifyDict[aId] = 1
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb(), self.getColumn()))
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)