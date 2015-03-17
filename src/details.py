import re
from itertools import izip
from collections import defaultdict, Counter

from jobTree.src.bioio import logger, reverseComplement

from lib.general_lib import formatRatio
from src.abstractClassifier import AbstractClassifier

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
from src.indels import deletionIterator, insertionIterator, firstIndel, codonPairIterator


class CodingInsertions(AbstractClassifier):
    """

    Does the alignment introduce insertions to the target genome?

    Reports a BED record for each such insertion

    Target insertion:

    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAnnotationDict()
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            insertions = list(insertionIterator(a, t, aln, mult3))
            if len(insertions) > 0:
                valueDict[aId] = [seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()) for start, stop, size in insertions if t.chromosomeCoordinateToCds(start) != None and t.chromosomeCoordinateToCds(stop) != None]
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class CodingMult3Insertions(CodingInsertions):
    """

    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.

    """

    def run(self):
        CodingInsertions.run(self, mult3=True)


class CodingDeletions(AbstractClassifier):
    """

    Does the alignment introduce deletions that are not a multiple of 3 to the target genome?

    Reports a BED record for each deletion. Since deletions are obviously not visible
    on the target genome, just has a 1 base record at this position.

    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA
             012345  67
    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            deletions = list(deletionIterator(a, t, aln, mult3))
            if len(deletions) > 0:
                valueDict[aId] = [seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()) for start, stop, size in deletions if t.chromosomeCoordinateToCds(start) != None]
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class CodingMult3Deletions(CodingDeletions):
    """

    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.

    """

    def run(self):
        CodingDeletions.run(self, mult3=True)


class FrameShift(AbstractClassifier):
    """

    Frameshifts are caused by coding indels that are not a multiple of 3. Reports a BED entry
    for the first frameshifting indel in this alignment.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        self.getAnnotationDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            transcript = self.transcriptDict[aId]
            annotatedTranscript = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            f = firstIndel(annotatedTranscript, transcript, aln, mult3=False, inversion=False)
            if f is not None:
                start, stop, size = f
                valueDict[aId] = seq_lib.chromosomeCoordinateToBed(transcript, start, stop, self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class AlignmentAbutsLeft(AbstractClassifier):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    If so, reports BED of entire transcript

    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
            elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class AlignmentAbutsRight(AbstractClassifier):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    If so, reports BED of entire transcript

    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
            elif aln.strand == "-" and aln.tStart == 0 and aln.qStart != 0:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class AlignmentPartialMap(AbstractClassifier):
    """

    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    If so, reports the entire transcript

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.qSize != aln.qEnd - aln.qStart:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class BadFrame(AbstractClassifier):
    """

    Looks for CDS sequences that are not a multiple of 3

    Will report a BED record of the transcript if true

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.getCdsLength() % 3 != 0:
                valueDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class BeginStart(AbstractClassifier):
    """

    Does the annotated CDS have a start codon (ATG) in the first 3 bases?

    Returns a BED record of the first 3 bases if this is NOT true

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.thickStart == t.thickStop == 0 or t.thickStop - t.thickStart < 3:
                continue
            s = t.getCds(self.seqDict)
            if not s.startswith("ATG"):
                valueDict[aId] = seq_lib.transcriptCoordinateToBed(t, 0, 2, self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class CdsGap(AbstractClassifier):
    """

    Are any of the CDS introns too short? Too short default is 30 bases.

    Reports a BED record for each intron interval that is too short.

    If mult3 is true, will only report on multiple of 3 gaps.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["alignment"]

    def mult3(self, t, shortIntronSize):
        records = []
        # only report if CdsGap is a multiple of 3
        if t.chromosomeInterval.strand is False:
            intronIntervals = t.intronIntervals[::-1]
        else:
            intronIntervals = t.intronIntervals
        for i in xrange(len(intronIntervals)):
            #is this intron coding?
            if t.exons[i].containsCds() is True and t.exons[i + 1].containsCds() is True:
                if len(intronIntervals[i]) <= shortIntronSize:
                    if len(intronIntervals[i]) % 3 == 0:
                        records.append(seq_lib.intervalToBed(t, intronIntervals[i], self.rgb(), self.getColumn()))
        if len(records) > 0:
            return records

    def notMult3(self, t, shortIntronSize):
        records = []
        if t.chromosomeInterval.strand is False:
            intronIntervals = t.intronIntervals[::-1]
        else:
            intronIntervals = t.intronIntervals
        # only report if CdsGap is a multiple of 3
        for i in xrange(len(intronIntervals)):
            #is this intron coding?
            if t.exons[i].containsCds() is True and t.exons[i + 1].containsCds() is True:
                if len(intronIntervals[i]) <= shortIntronSize:
                    if len(intronIntervals[i]) % 3 != 0:
                        records.append(seq_lib.intervalToBed(t, intronIntervals[i], self.rgb(), self.getColumn()))
        if len(records) > 0:
            return records

    def run(self, mult3=False, shortIntronSize=30):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = {}
        r = None
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            if mult3 is True:
                r = self.mult3(self.transcriptDict[aId], shortIntronSize)
            else:
                r = self.notMult3(self.transcriptDict[aId], shortIntronSize)
            if r is not None:
                valueDict[aId] = r
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class CdsMult3Gap(CdsGap):
    """

    See CdsGap for details. Runs it in mult3 mode.

    """

    def run(self, mult3=True, shortIntronSize=30):
        CdsGap.run(self, mult3, shortIntronSize)


class CdsNonCanonSplice(AbstractClassifier):
    """

    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    Reports two BED records of the four offending bases.

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["mutation"]

    def badSplice(self, donor, acceptor):
        m = {"GT": "AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def makeBed(self, t, intronInterval):
        tmp = []
        start, stop = intronInterval.start, intronInterval.start + 2
        tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
        start, stop = intronInterval.stop - 2, intronInterval.stop
        tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
        return tmp

    def run(self, shortIntronSize=30):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            for i, seq in enumerate(t.intronSequenceIterator(self.seqDict)):
                # make sure this intron is between coding exons
                if t.exons[i].containsCds() and t.exons[i + 1].containsCds():
                    if t.chromosomeInterval.strand is False:
                        if self.badSplice(reverseComplement(seq[:2]), reverseComplement(seq[-2:])) == True:
                            valueDict[aId] = self.makeBed(t, t.intronIntervals[i])
                    else:
                        if self.badSplice(seq[:2], seq[-2:]) == True:
                            valueDict[aId] = self.makeBed(t, t.intronIntervals[i])
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class CdsUnknownSplice(CdsNonCanonSplice):
    """

    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    subclasses cdsNonCanonSplice and just replaces the badSplice function

    reports 1 if TRUE, 0 if FALSE

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """

    def badSplice(self, donor, acceptor):
        m = {"GT": "AG", "GC": "AG", "AT": "AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        CdsNonCanonSplice.run(self)


class EndStop(AbstractClassifier):
    """

    Looks at the end of the coding region (thickEnd) and sees if the last
    three bases are a stop codon ('TAA', 'TGA', 'TAG')

    If this is NOT true, will report a BED record of the last 3 bases.
    
    Value will be NULL if there is insufficient information, which is defined as:
        1) thickStop - thickStart <= 9: (no useful CDS annotation)
        2) this alignment was not trans-mapped

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        stopCodons = ('TAA', 'TGA', 'TAG')
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.getCdsLength() <= 9:
                continue
            cds = t.getCds(self.seqDict)
            if cds[-3:].upper() not in stopCodons:
                valueDict[aId] = seq_lib.cdsCoordinateToBed(t, len(cds) - 3, len(cds), self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class InFrameStop(AbstractClassifier):
    """

    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 3 codons.

    Returns a BED record of the position of the in frame stop if it exists.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            # make sure this transcript has CDS
            #and more than 2 codons - can't have in frame stop without that
            if t.getCdsLength() >= 9:
                for i in xrange(9, t.getCdsLength() - 3, 3):
                    c = t.cdsCoordinateToAminoAcid(i, self.seqDict)
                    if c == "*":
                        valueDict[aId] = seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class NoCds(AbstractClassifier):
    """

    Looks to see if this transcript actually has a CDS, which is defined as having a
    thickStop-thickStart region of at least 1 codon. Adjusting cdsCutoff can change this.

    If True, reports entire transcript.

    Only reports if the original transcript had a CDS.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["alignment"]

    def run(self, cdsCutoff=3):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.getCdsLength() < cdsCutoff:
                a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
                if a.getCdsLength() > cdsCutoff:
                    valueDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class MinimumCdsSize(NoCds):
    """

    The smallest ORFs in any species are >10AA. So, we will flag any CDS smaller than this.

    Inherits NoCds and modifies cdsCutoff to do this.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        NoCds.run(self, cdsCutoff=30)


class ScaffoldGap(AbstractClassifier):
    """

    Does this alignment span a scaffold gap? (Defined as a 100bp run of Ns)

    Reports the entire alignment if it spans a scaffold gap

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getSeqDict()
        self.getTranscriptDict()
        valueDict = {}
        r = re.compile("[atgcATGC][N]{100}[atgcATGC]")
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            destSeq = self.seqDict[aln.tName][aln.tStart: aln.tEnd].upper()
            if re.search(r, destSeq) is not None:
                t = self.transcriptDict[aId]
                valueDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class UnknownBases(AbstractClassifier):
    """

    Does this alignment contain Ns in the target genome?

    Only looks mRNA bases, and restricts to CDS if cds is True

    Reports the entire transcript if True, regardless of cds flag

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["assembly"]

    def run(self, cds=False):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        r = re.compile("N+")
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if cds is True:
                s = t.getCds(self.seqDict)
                for m in re.finditer(r, s):
                        valueDict[aId] = seq_lib.cdsCoordinateToBed(t, m.start(), m.end(), self.rgb(), self.getColumn())
            else:
                s = t.getMRna(self.seqDict)
                for m in re.finditer(r, s):
                        valueDict[aId] = seq_lib.transcriptCoordinateToBed(t, m.start(), m.end(), self.rgb(), self.getColumn())
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class UnknownCdsBases(UnknownBases):
    """

    Inherits Unknown Bases and sets the cds flag to True.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def run(self):
        UnknownBases.run(self, cds=True)


class UtrGap(AbstractClassifier):
    """

    Are any UTR introns too short? Too short is defined as less than 30bp

    Reports on all such introns.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.chromosomeInterval.strand is False:
                intronIntervals = t.intronIntervals[::-1]
            else:
                intronIntervals = t.intronIntervals
            for i in xrange(len(intronIntervals)):
                if t.exons[i].containsCds() is False and t.exons[i + 1].containsCds() is False:
                    if len(intronIntervals[i]) <= shortIntronSize:
                        valueDict[aId].append(
                            seq_lib.intervalToBed(t, intronIntervals[i], self.rgb(), self.getColumn()))
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class UtrNonCanonSplice(AbstractClassifier):
    """

    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    Reports BED records for the bases that are bad.

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    TODO: this class is nearly identical to CdsNonCanonSplice. Devise a way to merge.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def rgb(self):
        return self.colors["generic"]

    def badSplice(self, donor, acceptor):
        m = {"GT": "AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def makeBed(self, t, intronInterval):
        tmp = []
        start, stop = intronInterval.start, intronInterval.start + 2
        tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
        start, stop = intronInterval.stop - 2, intronInterval.stop
        tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
        return tmp

    def run(self, shortIntronSize=30):
        logger.info("Starting detailed analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            for i, seq in enumerate(t.intronSequenceIterator(self.seqDict)):
                if not t.exons[i].containsCds() and not t.exons[i + 1].containsCds():
                    if t.chromosomeInterval.strand is False:
                        if self.badSplice(reverseComplement(seq[:2]), reverseComplement(seq[-2:])) == True:
                            valueDict[aId] = self.makeBed(t, t.intronIntervals[i])
                    else:
                        if self.badSplice(seq[:2], seq[-2:]) == True:
                            valueDict[aId] = self.makeBed(t, t.intronIntervals[i])
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class UtrUnknownSplice(UtrNonCanonSplice):
    """

    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    subclasses CdsNonCanonSplice and just replaces the badSplice function

    reports 1 if TRUE, 0 if FALSE

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """

    @staticmethod
    def _getType():
        return "TEXT"

    def badSplice(self, donor, acceptor):
        m = {"GT": "AG", "GC": "AG", "AT": "AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        UtrNonCanonSplice.run(self)


class Nonsynonymous(AbstractClassifier):
    """

    Do any base changes introduce nonsynonmous changes?

    This is filtered for codons within frameshifts, unless another frameshift restores the frame.

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def rgb(self):
        return self.colors["nonsynon"]

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        self.getRefTwoBit()
        self.getAlignmentDict()
        valueDict = defaultdict(list)
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refTwoBit):
                if target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) != seq_lib.codonToAminoAcid(query_codon):
                    valueDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb(), self.getColumn()))
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class Synonymous(AbstractClassifier):
    """

    Do any base changes introduce nonsynonmous changes?

    This is filtered for codons within frameshifts, unless another frameshift restores the frame.

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def rgb(self):
        return self.colors["synon"]        

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getAnnotationDict()
        self.getSeqDict()
        self.getRefTwoBit()
        self.getAlignmentDict()
        valueDict = defaultdict(list)
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            a = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refTwoBit):
                if target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) == seq_lib.codonToAminoAcid(query_codon):
                        valueDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb(), self.getColumn()))
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)


class Paralogy(AbstractClassifier):
    """

    Does this transcript appear more than once in the transcript dict?

    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def rgb(self):
        return self.colors["mutation"]       

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        counts = Counter(psl_lib.removeAlignmentNumber(aId) for aId in self.aIds if aId in self.transcriptDict)
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            elif counts[psl_lib.removeAlignmentNumber(aId)] > 1:
                t = self.transcriptDict[aId]
                valueDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn() +"_{}_Copies".format(counts[psl_lib.removeAlignmentNumber(aId)] - 1))
        logger.info(
            "Details {} on {} is finished. {} records failed".format(self.genome, self.getColumn(), len(valueDict)))
        self.dumpValueDict(valueDict)        
