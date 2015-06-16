import re
from itertools import izip
from collections import defaultdict, Counter

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

from src.abstractClassifier import AbstractAugustusClassifier


def gapMergeExons(t, shortIntronSize=30):
    """
    takes a transcript and merges all introns shorter than shortIntronSize, returning a list of intervals
    that represent the exon boundaries
    """


class AugustusParalogy(AbstractAugustusClassifier):
    """
    Does this transcript appear more than once in the augustus transcript dict?
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        self.getAugustusTranscriptDict()
        counts = Counter(psl_lib.removeAugustusAlignmentNumber(aug_aId) for aug_aId in self.augustusTranscriptDict)
        detailsDict = {}
        classifyDict = {}
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            if counts[psl_lib.removeAugustusAlignmentNumber(aug_aId)] > 1:
                detailsDict[aug_aId] = seq_lib.transcriptToBed(aug_t, self.rgb, self.column + "_{}_Copies".format(
                    counts[psl_lib.removeAugustusAlignmentNumber(aug_aId)] - 1))
                classifyDict[aug_aId] = 1
            else:
                classifyDict[aug_aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AugustusExonGain(AbstractAugustusClassifier):
    """
    Does the augustus version of this transcript add an exon?
    This is calculated by looking at the exon boundary intervals between the genePreds
    we expect new exons to not be within the gap-merged exons of the transMap genes.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        self.getAugustusTranscriptDict()
        self.getTranscriptDict()
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            t = self.transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
