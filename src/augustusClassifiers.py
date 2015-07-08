from collections import defaultdict, Counter
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from src.abstractClassifier import AbstractAugustusClassifier


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
                                                               counts[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
                                                               - 1))
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

    def run(self, shortIntronSize=30):
        self.getAugustusTranscriptDict()
        self.getTranscriptDict()
        classifyDict = {}
        detailsDict = defaultdict(list)
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            if psl_lib.removeAugustusAlignmentNumber(aug_aId) not in self.transcriptDict:
                continue
            t = self.transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
            if aug_t.strand != t.strand or aug_t.chromosome != t.chromosome:
                continue
            aug_t_intervals = aug_t.exonIntervals
            merged_t_intervals = seq_lib.gapMergeIntervals(t.exonIntervals, gap=shortIntronSize)
            for interval in aug_t_intervals:
                if seq_lib.intervalNotIntersectIntervals(merged_t_intervals, interval):
                    classifyDict[aug_aId] = 1
                    detailsDict[aug_aId].append(interval.getBed(self.rgb, self.column))
            if aug_aId not in classifyDict:
                classifyDict[aug_aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AugustusExonLoss(AbstractAugustusClassifier):
    """
    Does the augustus version of this transcript lose an exon?
    This is calculated by looking at the exon boundary intervals between the genePreds
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30):
        self.getAugustusTranscriptDict()
        self.getTranscriptDict()
        classifyDict = {}
        detailsDict = defaultdict(list)
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            if psl_lib.removeAugustusAlignmentNumber(aug_aId) not in self.transcriptDict:
                continue
            t = self.transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
            if aug_t.strand != t.strand or aug_t.chromosome != t.chromosome:
                continue
            aug_t_intervals = aug_t.exonIntervals
            merged_t_intervals = seq_lib.gapMergeIntervals(t.exonIntervals, gap=shortIntronSize)
            for interval in merged_t_intervals:
                if seq_lib.intervalNotIntersectIntervals(aug_t_intervals, interval):
                    classifyDict[aug_aId] = 1
                    detailsDict[aug_aId].append(interval.getBed(self.rgb, self.column))
            if aug_aId not in classifyDict:
                classifyDict[aug_aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AugustusNotSimilarExonBoundaries(AbstractAugustusClassifier):
    """
    Does the augustus version of this transcript have the same exon boundaries, within a wiggle room range?
    Returns True if any exons are NOT in the same boundaries, and also reports each such splice junction
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30, wiggleRoom=10):
        self.getAugustusTranscriptDict()
        self.getTranscriptDict()
        classifyDict = {}
        detailsDict = defaultdict(list)
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            if psl_lib.removeAugustusAlignmentNumber(aug_aId) not in self.transcriptDict:
                continue
            t = self.transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
            if aug_t.strand != t.strand or aug_t.chromosome != t.chromosome:
                continue
            merged_t_intervals = seq_lib.gapMergeIntervals(t.exonIntervals, gap=shortIntronSize)
            aug_t_intervals = aug_t.exonIntervals
            for interval in merged_t_intervals:
                if seq_lib.intervalNotWithinWiggleRoomIntervals(aug_t_intervals, interval, wiggleRoom):
                    classifyDict[aug_aId] = 1
                    detailsDict[aug_aId].append(interval.getBed(self.rgb, self.column))
            if aug_aId not in classifyDict:
                classifyDict[aug_aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AugustusSameStartStop(AbstractAugustusClassifier):
    """
    Does the augustus transcript have the exact same start bases as the transMap transcript?
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        self.getAugustusTranscriptDict()
        self.getTranscriptDict()
        classifyDict = {}
        detailsDict = {}
        for aug_aId, aug_t in self.augustusTranscriptDict.iteritems():
            if psl_lib.removeAugustusAlignmentNumber(aug_aId) not in self.transcriptDict:
                continue
            t = self.transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
            if aug_t.strand != t.strand or aug_t.chromosome != t.chromosome or t.thickStart == t.thickStop:
                continue
            if t.thickStart == aug_t.thickStart and t.thickStop == aug_t.thickStop:
                classifyDict[aug_aId] = 0
                s = t.getCdsLength()
                detailsDict[aug_aId] = [seq_lib.cdsCoordinateToBed(t, 0, 3, self.rgb, self.column),
                                        seq_lib.cdsCoordinateToBed(t, s - 3, s, self.rgb, self.column)]
            else:
                classifyDict[aug_aId] = 1
        self.dumpValueDicts(classifyDict, detailsDict)
