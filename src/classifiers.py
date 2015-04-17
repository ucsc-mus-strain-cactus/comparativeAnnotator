import re
from itertools import izip
from collections import defaultdict, Counter

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

from jobTree.src.bioio import logger, reverseComplement

from src.abstractClassifier import AbstractClassifier
from src.helperFunctions import deletionIterator, insertionIterator, frameShiftIterator, codonPairIterator, compareIntronToReference

class CodingInsertions(AbstractClassifier):
    """
    Does the alignment introduce insertions to the target genome?

    Reports a BED record for each such insertion

    Target insertion:

    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            insertions = [seq_lib.chromosomeRegionToBed(t, start, stop, self.rgb, self.column) for start, stop, size in insertionIterator(a, t, aln, mult3) if start >= t.thickStart and stop < t.thickStop]
            if len(insertions) > 0:
                detailsDict[aId] = insertions
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CodingMult3Insertions(CodingInsertions):
    """
    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingInsertions.run(self, mult3=True)


class CodingDeletions(AbstractClassifier):
    """
    Does the alignment introduce deletions that are not a multiple of 3 to the target genome?

    Reports a BED record for each deletion. Since deletions are not visible on the target genome, just has a 0 base 
    record at this position.

    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            deletions = [seq_lib.chromosomeRegionToBed(t, start, stop, self.rgb, self.column) for start, stop, size in deletionIterator(a, t, aln, mult3) if start >= t.thickStart and stop < t.thickStop]            
            if len(deletions) > 0:
                detailsDict[aId] = deletions
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CodingMult3Deletions(CodingDeletions):
    """
    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingDeletions.run(self, mult3=True)


class StartOutOfFrame(AbstractClassifier):
    """
    StartOutOfFrame are caused when the starting CDS base of the lifted over transcript is not in the original frame.
    If True, reports the first 3 bases.

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            aln = self.alignmentDict[aId]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            if a.transcriptCoordinateToCds(aln.targetCoordinateToQuery(t.cdsCoordinateToChromosome(0))) % 3 != 0:
                classifyDict[aId] = 1
                detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, 0, 3, self.rgb, self.column)
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class FrameShift(AbstractClassifier):
    """
    Frameshifts are caused by coding indels that are not a multiple of 3. Reports a BED entry
    spanning all blocks of coding bases that are frame-shifted. Must have at least 25 codons.

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            # do not include noncoding transcripts or lift-overs that contain less than 1 codon
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            frame_shifts = list(frameShiftIterator(a, t, aln))
            if len(frame_shifts) == 0:
                classifyDict[aId] = 0
                continue
            indel_starts, indel_stops, spans = zip(*frame_shifts)
            # calculate cumulative frame by adding each span and taking mod 3 - zeroes imply regaining frame
            # note that this code prepends a 0 to the list, offsetting all values by 1. This is useful.
            cumulative_frame = map(lambda x: x % 3, reduce(lambda l, v: (l.append(l[-1] + v) or l), spans, [0]))
            # every start is when a zero existed in the previous spot in cumulative_frame
            windowed_starts = [x for x, y in izip(indel_starts, cumulative_frame) if y == 0 or x == indel_starts[0]]
            # every stop is when a zero exists at this cumulative_frame
            windowed_stops = [x for x, y in izip(indel_stops, cumulative_frame[1:]) if y == 0]
            # sanity check
            assert any([len(windowed_starts) == len(windowed_stops), len(windowed_starts) - 1 == len(windowed_stops)]), (self.genome, self.column, aId)
            # now we need to fix frame and stops - if this shift extends to the end of the transcript, add that stop
            # additionally, if this is a negative strand transcript, flip starts/stops so that start is always < stop
            if len(windowed_stops) < len(windowed_starts) and t.strand is False:
                windowed_stops.append(t.thickStart)
                windowed_stops, windowed_starts = windowed_starts, windowed_stops
            elif len(windowed_stops) < len(windowed_starts):
                windowed_stops.append(t.thickStop)
            elif t.strand is False:
                windowed_stops, windowed_starts = windowed_starts, windowed_stops
            for start, stop in izip(windowed_starts, windowed_stops):
                detailsDict[aId].append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb, self.column))
            classifyDict[aId] = 1
        self.dumpValueDicts(classifyDict, detailsDict)


class AlignmentAbutsLeft(AbstractClassifier):
    """
    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    If so, reports BED of entire transcript

    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....

    Doesn't need a check for pre-existing because that doesn't matter.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb, self.column)
                classifyDict[aId] = 1
            elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AlignmentAbutsRight(AbstractClassifier):
    """
    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    If so, reports BED of entire transcript

    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|

    Doesn't need a check for pre-existing because that doesn't matter.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb, self.column)
                classifyDict[aId] = 1
            elif aln.strand == "-" and aln.tStart == 0 and aln.qStart != 0:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AlignmentPartialMap(AbstractClassifier):
    """
    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    If so, reports the entire transcript

    Doesn't need a check for pre-existing because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.qSize != aln.qEnd - aln.qStart:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class BadFrame(AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least 25 codons.

    Will report a BED record of the transcript if true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            if t.getCdsLength() % 3 != 0 and a.getCdsLength() % 3 != 0:
                detailsDict[aId] = seq_lib.chromosomeCoordinateToBed(t, t.thickStart, t.thickStop, self.colors["input"], self.column)
                classifyDict[aId] = 1
            elif t.getCdsLength() % 3 != 0:
                detailsDict[aId] = seq_lib.chromosomeCoordinateToBed(t, t.thickStart, t.thickStop, self.rgb, self.column)
                classifyDict[aId] = 1                
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class BeginStart(AbstractClassifier):
    """
    Does the lifted over CDS have the same 3 start bases as the original transcript?
    AND are these bases 'ATG'? Must have at least 25 codons.

    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getAlignmentDict()
        self.getAnnotationDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            aln = self.alignmentDict[aId]
            # do not include noncoding transcripts or lift-overs that contain less than 1 codon
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            cds_positions = [t.chromosomeCoordinateToCds(aln.queryCoordinateToTarget(a.cdsCoordinateToTranscript(i))) for i in xrange(3)]
            if None in cds_positions:
                detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, 0, 3, self.rgb, self.column)
                classifyDict[aId] = 1
            elif t.getCds(self.seqDict)[:3] != "ATG":
                if a.getCds(self.refTwoBit)[:3] != "ATG":
                    detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, 0, 3, self.colors["input"], self.column)
                else:
                    detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, 0, 3, self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CdsGap(AbstractClassifier):
    """
    Are any of the CDS introns too short? Too short default is 30 bases.

    Reports a BED record for each intron interval that is too short.

    Not currently implementing checking for pre-existing because its so rare (none in VM4)
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for i, intron in enumerate(t.intronIntervals):
                if len(intron) >= shortIntronSize:
                    continue
                elif "N" in intron.getSequence(self.seqDict):
                    continue
                elif not (intron.start >= t.thickStart and intron.stop < t.thickStop):
                    continue
                detailsDict[aId].append(seq_lib.intervalToBed(t, intron, self.rgb, self.column))
                classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CdsMult3Gap(AbstractClassifier):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]    

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for i, intron in enumerate(t.intronIntervals):
                if len(intron) >= shortIntronSize:
                    continue
                elif len(intron) % 3 != 0:
                    continue
                elif "N" in intron.getSequence(self.seqDict):
                    continue
                elif not (intron.start >= t.thickStart and intron.stop < t.thickStop):
                    continue
                detailsDict[aId].append(seq_lib.intervalToBed(t, intron, self.rgb, self.column))
                classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UtrGap(AbstractClassifier):
    """
    Are any UTR introns too short? Too short is defined as less than 30bp

    Reports on all such introns.

    Not currently implementing checking for pre-existing because its so rare (397 in VM4)
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for i, intron in enumerate(t.intronIntervals):
                if len(intron) >= shortIntronSize:
                    continue
                elif "N" in intron.getSequence(self.seqDict):
                    continue
                elif intron.start >= t.thickStart and intron.stop < t.thickStop:
                    continue
                detailsDict[aId].append(seq_lib.intervalToBed(t, intron, self.rgb, self.column))
                classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UnknownGap(AbstractClassifier):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.

    Not implementing looking for pre-existing because that doesn't make sense.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for i, intron in enumerate(t.intronIntervals):
                if len(intron) >= shortIntronSize:
                    continue
                elif "N" not in intron.getSequence(self.seqDict):
                    continue
                detailsDict[aId].append(seq_lib.intervalToBed(t, intron, self.rgb, self.column))
                classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)        


class CdsNonCanonSplice(AbstractClassifier):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    Reports two BED records of the four offending bases.

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    canonical = {"GT": "AG"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for intron in t.intronIntervals:
                if len(intron) <= shortIntronSize:
                    continue
                elif not (intron.start >= t.thickStart and intron.stop < t.thickStop):
                    continue
                seq = intron.getSequence(self.seqDict, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.canonical or self.canonical[donor] != acceptor:
                    classifyDict[aId] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compareIntronToReference(intron, a, t, aln, self.canonical, self.refDict) is True:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.colors["input"], self.column))
                    else:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb, self.column))
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CdsUnknownSplice(AbstractClassifier):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for intron in t.intronIntervals:
                if len(intron) <= shortIntronSize:
                    continue
                elif not (intron.start >= t.thickStart and intron.stop < t.thickStop):
                    continue
                seq = intron.getSequence(self.seqDict, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.non_canonical or self.non_canonical[donor] != acceptor:
                    classifyDict[aId] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compareIntronToReference(intron, a, t, aln, self.non_canonical, self.refDict) is True:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.colors["input"], self.column))
                    else:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb, self.column))
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UtrNonCanonSplice(AbstractClassifier):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    canonical = {"GT": "AG"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for intron in t.intronIntervals:
                if len(intron) <= shortIntronSize:
                    continue
                elif intron.start >= t.thickStart and intron.stop < t.thickStop:
                    continue
                seq = intron.getSequence(self.seqDict, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.canonical or self.canonical[donor] != acceptor:
                    classifyDict[aId] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compareIntronToReference(intron, a, t, aln, self.canonical, self.refDict) is True:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.colors["input"], self.column))
                    else:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb, self.column))
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UtrUnknownSplice(AbstractClassifier):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            for intron in t.intronIntervals:
                if len(intron) <= shortIntronSize:
                    continue
                elif intron.start >= t.thickStart and intron.stop < t.thickStop:
                    continue
                seq = intron.getSequence(self.seqDict, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.non_canonical or self.non_canonical[donor] != acceptor:
                    classifyDict[aId] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compareIntronToReference(intron, a, t, aln, self.non_canonical, self.refDict) is True:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.colors["input"], self.column))
                    else:
                        detailsDict[aId].append(seq_lib.spliceIntronIntervalToBed(t, intron, self.rgb, self.column))
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class EndStop(AbstractClassifier):
    """
    Looks at the last three bases of the annotated transcript and determines if they properly match the last 3 bases 
    of the lifted over transcript AND that those bases are in ('TAA', 'TGA', 'TAG'). Must have at least 25 codons.

    If this is NOT true, will report a BED record of the last 3 bases.
    
    Value will be NULL if there is insufficient information, which is defined as:
        1) thickStop - thickStart <= 9: (no useful CDS annotation)
        2) this alignment was not trans-mapped
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        stopCodons = ('TAA', 'TGA', 'TAG')
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            aln = self.alignmentDict[aId]
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            s = t.getCdsLength()
            cds_positions = [t.chromosomeCoordinateToCds(aln.queryCoordinateToTarget(a.cdsCoordinateToTranscript(i))) for i in xrange(s - 4, s - 1)]
            if None in cds_positions or t.getCds(self.seqDict)[-3:] not in stopCodons:
                # does this problem exist in the reference?
                if a.getCds(self.refDict)[-3:] not in stopCodons:
                    detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, s - 3, s, self.colors["input"], self.column)
                else:
                    detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, s - 3, s, self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class InFrameStop(AbstractClassifier):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 25 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        self.getRefDict()
        self.getAlignmentDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            aln = self.alignmentDict[aId]
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            # TODO: this will miss an inframe stop if it is the last 3 bases that are not the annotated stop.
            # use the logic from EndStop to flag this
            codons = list(codonPairIterator(a, t, aln, self.seqDict, self.refDict))[:-1]
            for i, target_codon, query_codon in codons:
                if seq_lib.codonToAminoAcid(target_codon) == "*":
                    if target_codon == query_codon:
                        detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.colors["input"], self.column))
                    else:
                        detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb, self.column))
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class ShortCds(AbstractClassifier):
    """
    Looks to see if this transcript actually has a CDS, which is defined as having a CDS region of 
    at least 25 codons. Adjusting cdsCutoff can change this.

    If True, reports entire transcript.

    Only reports if the original transcript had a CDS.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, cdsCutoff=75):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            # do not include noncoding transcripts
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            if a.getCdsLength() < 3:
                continue
            elif a.getCdsLength() <= cdsCutoff:
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.colors["input"], self.column)
                classifyDict[aId] = 1                
            elif t.getCdsLength() <= cdsCutoff:
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class ScaffoldGap(AbstractClassifier):
    """
    Does this alignment span a scaffold gap? (Defined as a 100bp run of Ns)

    Reports the entire alignment if it spans a scaffold gap
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getAlignmentDict()
        self.getSeqDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        r = re.compile("[atgcATGC][nN]{100}[atgcATGC]")
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            destSeq = self.seqDict[aln.tName][aln.tStart:aln.tEnd]
            if re.search(r, destSeq) is not None:
                t = self.transcriptDict[aId]
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb, self.column)
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UnknownBases(AbstractClassifier):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, cds=False):
        logger.info("Starting analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        r = re.compile("N+")
        for aId, t in self.transcriptDict.iteritems():
            if cds is True:
                s = t.getCds(self.seqDict)
                tmp = [seq_lib.cdsCoordinateToBed(t, m.start(), m.end(), self.rgb, self.column) for m in re.finditer(r, s)]
            else:
                s = t.getMRna(self.seqDict)
                tmp = [seq_lib.transcriptCoordinateToBed(t, m.start(), m.end(), self.rgb, self.column) for m in re.finditer(r, s)]
            if len(tmp) > 0:
                detailsDict[aId] = tmp
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UnknownCdsBases(UnknownBases):
    """
    Inherits Unknown Bases and sets the cds flag to True.
    """

    def run(self):
        UnknownBases.run(self, cds=True)


class Nonsynonymous(AbstractClassifier):
    """
    Do any base changes introduce nonsynonmous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["nonsynon"]

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        self.getRefDict()
        self.getAlignmentDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refDict):
                if "N" not in target_codon and target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) != seq_lib.codonToAminoAcid(query_codon):
                    detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb, self.column))
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class Synonymous(AbstractClassifier):
    """
    Do any base changes introduce nonsynonmous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["synon"]        

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        self.getRefDict()
        self.getAlignmentDict()
        detailsDict = defaultdict(list)
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            if a.getCdsLength() <= 75 or t.getCdsLength() <= 75:
                continue
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refDict):
                if target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) == seq_lib.codonToAminoAcid(query_codon):
                    detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb, self.column))
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class Paralogy(AbstractClassifier):
    """
    Does this transcript appear more than once in the transcript dict?
    """
    @property
    def rgb(self):
        return self.colors["mutation"]       

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.column, self.genome))
        self.getTranscriptDict()
        counts = Counter(psl_lib.removeAlignmentNumber(aId) for aId in self.transcriptDict)
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if counts[psl_lib.removeAlignmentNumber(aId)] > 1:
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb, self.column + "_{}_Copies".format(counts[psl_lib.removeAlignmentNumber(aId)] - 1))
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)
