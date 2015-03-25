import re
from itertools import izip
from collections import defaultdict, Counter

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

from jobTree.src.bioio import logger, reverseComplement

from src.abstractClassifier import AbstractClassifier
from src.helperClasses import GapFinder, SpliceSiteAnalysis
from src.helperFunctions import deletionIterator, insertionIterator, frameShiftIterator, codonPairIterator
from src.helperFunctions import rearrangementIterator

class CodingInsertions(AbstractClassifier):
    """
    Does the alignment introduce insertions to the target genome?

    Reports a BED record for each such insertion

    Target insertion:

    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
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
            insertions = [seq_lib.chromosomeRegionToBed(t, start, stop, self.rgb(), self.getColumn()) for start, stop, \
                          size in insertionIterator(a, t, aln, mult3, rearrangement=False) if start >= t.thickStart \
                           and stop <= t.thickStop]
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
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
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
            deletions = [seq_lib.chromosomeRegionToBed(t, start, stop, self.rgb(), self.getColumn()) for start, stop, 
                         size in deletionIterator(a, t, aln, mult3, rearrangement=False) if start >= t.thickStart and stop \
                         <= t.thickStop]            
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


class Rearrangements(AbstractClassifier):
    """

    Are there any rearrangements in these alignments? Will show both insertion and deletion style rearrangements

    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
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
            rearrangements = [seq_lib.chromosomeRegionToBed(t, start, stop, self.rgb(), self.getColumn()) for start, stop,
                          size in rearrangementIterator(a, t, aln)]
            if len(rearrangements) > 0:
                detailsDict[aId] = rearrangements
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class FrameMismatch(AbstractClassifier):
    """
    Frameshifts are caused by coding indels that are not a multiple of 3. Reports a BED entry
    spanning all blocks of coding bases that are frame-shifted.
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
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
            s = list(frameShiftIterator(a, t, aln))
            if len(s) == 0:
                classifyDict[aId] = 0
                continue
            elif len(s) == 1:
                start, stop, size = s[0]
                detailsDict[aId] = seq_lib.chromosomeCoordinateToBed(t, start, t.stop, self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                tmp = []
                for i in xrange(1, len(s), 2):
                    start = s[i-1][0]
                    stop = s[i][1]
                    tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
                if i % 2 == 0 and i < len(s):
                    start = s[-1][0]
                    tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, t.thickStop, self.rgb(), self.getColumn()))
                detailsDict[aId] = tmp
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
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
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
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.strand == "+" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            elif aln.strand == "-" and aln.tStart == 0 and aln.qStart != 0:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class AlignmentPartialMap(AbstractClassifier):
    """
    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    If so, reports the entire transcript
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            if aln.qSize != aln.qEnd - aln.qStart:
                detailsDict[aId] = seq_lib.transcriptToBed(self.transcriptDict[aId], self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class BadFrame(AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3

    Will report a BED record of the transcript if true
    """
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getAlignmentDict()
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if t.getCdsLength() % 3 != 0:
                detailsDict[aId] = seq_lib.chromosomeCoordinateToBed(t, t.thickStart, t.thickStop, self.rgb(),
                                                                     self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class BeginStart(AbstractClassifier):
    """
    Does the annotated CDS have a start codon (ATG) in the first 3 bases?

    Returns a BED record of the first 3 bases if this is NOT true
    """
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if t.thickStart == t.thickStop == 0 or t.thickStop - t.thickStart < 3:
                continue
            s = t.getCds(self.seqDict)
            if not s.startswith("ATG"):
                detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, 0, 2, self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class CdsGap(GapFinder):
    """
    Are any of the CDS introns too short? Too short default is 30 bases.

    Reports a BED record for each intron interval that is too short.

    If mult3 is true, will only report on multiple of 3 gaps.
    """
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30, mult3=False, coding=True):
        GapFinder.run(self, shortIntronSize=shortIntronSize, mult3=mult3, coding=coding)


class CdsMult3Gap(GapFinder):
    """
    See CdsGap for details. Runs it in mult3 mode.
    """
    def rgb(self):
        return self.colors["alignment"]    

    def run(self, shortIntronSize=30, mult3=True, coding=True):
        GapFinder.run(self, shortIntronSize=shortIntronSize, mult3=mult3, coding=coding)


class UtrGap(GapFinder):
    """
    Are any UTR introns too short? Too short is defined as less than 30bp

    Reports on all such introns.
    """
    def rgb(self):
        return self.colors["alignment"]

    def run(self, shortIntronSize=30, mult3=False, coding=False):
        GapFinder.run(self, shortIntronSize=shortIntronSize, mult3=mult3, coding=coding)


class CdsNonCanonSplice(SpliceSiteAnalysis):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    Reports two BED records of the four offending bases.

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30, coding=True, canonical=True):
        SpliceSiteAnalysis.run(self, canonical=canonical, coding=coding, shortIntronSize=shortIntronSize)


class CdsUnknownSplice(SpliceSiteAnalysis):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self, shortIntronSize=30, coding=True, canonical=False):
        SpliceSiteAnalysis.run(self, canonical=canonical, coding=coding, shortIntronSize=shortIntronSize)


class UtrNonCanonSplice(SpliceSiteAnalysis):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self, shortIntronSize=30, coding=False, canonical=True):
        SpliceSiteAnalysis.run(self, canonical=canonical, coding=coding, shortIntronSize=shortIntronSize)


class UtrUnknownSplice(SpliceSiteAnalysis):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self, shortIntronSize=30, coding=False, canonical=False):
        SpliceSiteAnalysis.run(self, canonical=canonical, coding=coding, shortIntronSize=shortIntronSize)


class EndStop(AbstractClassifier):
    """
    Looks at the end of the coding region (thickEnd) and sees if the last
    three bases are a stop codon ('TAA', 'TGA', 'TAG')

    If this is NOT true, will report a BED record of the last 3 bases.
    
    Value will be NULL if there is insufficient information, which is defined as:
        1) thickStop - thickStart <= 9: (no useful CDS annotation)
        2) this alignment was not trans-mapped
    """
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        stopCodons = ('TAA', 'TGA', 'TAG')
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if t.getCdsLength() <= 9:
                continue
            cds = t.getCds(self.seqDict)
            if cds[-3:] not in stopCodons:
                detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, len(cds) - 3, len(cds), self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class InFrameStop(AbstractClassifier):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 3 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            # make sure this transcript has CDS
            #and more than 2 codons - can't have in frame stop without that
            if t.getCdsLength() < 9:
                continue
            for i in xrange(9, t.getCdsLength() - 3, 3):
                c = t.cdsCoordinateToAminoAcid(i, self.seqDict)
                if c == "*":
                    detailsDict[aId] = seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb(), self.getColumn())
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class NoCds(AbstractClassifier):
    """
    Looks to see if this transcript actually has a CDS, which is defined as having a
    thickStop-thickStart region of at least 10 codons. Adjusting cdsCutoff can change this.

    If True, reports entire transcript.

    Only reports if the original transcript had a CDS.
    """
    def rgb(self):
        return self.colors["alignment"]

    def run(self, cdsCutoff=30):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if t.getCdsLength() < cdsCutoff:
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class ScaffoldGap(AbstractClassifier):
    """
    Does this alignment span a scaffold gap? (Defined as a 100bp run of Ns)

    Reports the entire alignment if it spans a scaffold gap
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
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
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn())
                classifyDict[aId] = 1
            else:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class UnknownBases(AbstractClassifier):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    def rgb(self):
        return self.colors["assembly"]

    def run(self, cds=False):
        logger.info("Starting analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        self.getSeqDict()
        detailsDict = {}
        classifyDict = {}
        r = re.compile("N+")
        for aId, t in self.transcriptDict.iteritems():
            if cds is True:
                s = t.getCds(self.seqDict)
                tmp = [seq_lib.cdsCoordinateToBed(t, m.start(), m.end(), self.rgb(), self.getColumn()) for m in \
                       re.finditer(r, s)]
            else:
                s = t.getMRna(self.seqDict)
                tmp = [seq_lib.transcriptCoordinateToBed(t, m.start(), m.end(), self.rgb(), self.getColumn()) for m in \
                       re.finditer(r, s)]
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
    Do any base changes introduce nonsynonmous changes?

    This is filtered for codons within frameshifts, unless another frameshift restores the frame.
    """
    def rgb(self):
        return self.colors["nonsynon"]

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
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
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refDict):
                if target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) != \
                                                                                  seq_lib.codonToAminoAcid(query_codon):
                    detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i - 3, i, self.rgb(), self.getColumn()))
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class Synonymous(AbstractClassifier):
    """
    Do any base changes introduce nonsynonmous changes?

    This is filtered for codons within frameshifts, unless another frameshift restores the frame.
    """
    def rgb(self):
        return self.colors["synon"]        

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
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
            for i, target_codon, query_codon in codonPairIterator(a, t, aln, self.seqDict, self.refDict):
                if target_codon != query_codon and seq_lib.codonToAminoAcid(target_codon) == \
                                                                                  seq_lib.codonToAminoAcid(query_codon):
                    detailsDict[aId].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, self.rgb(), self.getColumn()))
                    classifyDict[aId] = 1
            if aId not in classifyDict:
                classifyDict[aId] = 0
        self.dumpValueDicts(classifyDict, detailsDict)


class Paralogy(AbstractClassifier):
    """
    Does this transcript appear more than once in the transcript dict?
    """
    def rgb(self):
        return self.colors["mutation"]       

    def run(self):
        logger.info("Starting details analysis {} on {}".format(self.getColumn(), self.genome))
        self.getTranscriptDict()
        counts = Counter(psl_lib.removeAlignmentNumber(aId) for aId in self.aIds if aId in self.transcriptDict)
        detailsDict = {}
        classifyDict = {}
        for aId, t in self.transcriptDict.iteritems():
            if counts[psl_lib.removeAlignmentNumber(aId)] > 1:
                detailsDict[aId] = seq_lib.transcriptToBed(t, self.rgb(), self.getColumn() + "_{}_Copies".format( \
                                                           counts[psl_lib.removeAlignmentNumber(aId)] - 1))
        self.dumpValueDicts(classifyDict, detailsDict)        
