from src.abstractClassifier import AbstractClassifier

class GapFinder(AbstractClassifier):
    """
    Finds improperly sized UTRs. 
    """

    def gapFinder(self, t, shortIntronSize, mult3, coding):
        """
        Looks for introns that are smaller than shortIntronSize.
        Reports on these depending on the mult3 and coding flags
        mult3: True for only mult3, False for no mult3, None for either.
        coding: True for only coding, False for non-coding, None for either
        """
        records = []
        if t.chromosomeInterval.strand is False:
            intronIntervals = reversed(t.intronIntervals)
        else:
            intronIntervals = t.intronIntervals
        for i, intron in enumerate(intronIntervals):
            if len(intron) <= shortIntronSize:
                if intron.start >= t.thickStart and intron.stop <= t.thickStop and (coding is True or coding is None):
                    if len(intron) % 3 == 0 and (mult3 is True or mult3 is None):
                        records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.column()))
                    elif len(intron) % 3 != 0 and (mult3 is False or mult3 is None):
                        records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.column()))
                elif intron.start <= t.thickStart or intron.stop >= t.thickStop and (coding is False or coding is None):
                    if len(intron) % 3 == 0 and (mult3 is True or mult3 is None):
                        records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.column()))
                    elif len(intron) % 3 != 0 and (mult3 is False or mult3 is None):
                        records.append(seq_lib.intervalToBed(t, intron, self.rgb(), self.column()))
        return records

    def run(self, shortIntronSize=30, mult3=None, coding=None):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
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