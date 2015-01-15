import re
from itertools import izip

from lib.general_lib import formatRatio
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib


class AlignmentAbutsLeft(object):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
            return 1
        elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
            return 1
        return 0


class AlignmentAbutsRight(object):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        if aln.strand == "+" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
            return 1
        elif aln.strand == "-" and aln.tStart == 0 and aln.qStart != 0:
            return 1
        return 0


class AlignmentCoverage(object):
    """

    Calculates alignment coverage:

    (matches + mismatches) / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    """

    @staticmethod
    def _getType():
        return "REAL"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        r = formatRatio(aln.matches + aln.misMatches, aln.matches + aln.misMatches + aln.qNumInsert)
        assert r >= 0 and r <= 1
        return r


class AlignmentIdentity(object):
    """

    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    """

    @staticmethod
    def _getType():
        return "REAL"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        r = formatRatio(aln.matches, aln.matches + aln.misMatches + aln.qNumInsert)
        assert r >= 0 and r <= 1
        return r


class AlignmentPartialMap(object):
    """

    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        if aln.qSize != aln.qEnd - aln.qStart:
            return 1
        return 0


class BadFrame(object):
    """

    Looks for CDS sequences that are not a multiple of 3

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        if t.getCdsLength() % 3 != 0:
            return 1
        return 0


class TranscriptID(object):
    """
    Creates a column representing the transcript ID
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        return psl_lib.removeAlignmentNumber(aId)


class GeneID(object):
    """
    Creates a column representing the gene ID
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        return attributeDict[psl_lib.removeAlignmentNumber(aId)].geneID


class GeneName(object):
    """
    Creates a column representing the gene name
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        return attributeDict[psl_lib.removeAlignmentNumber(aId)].geneName


class GeneType(object):
    """
    Creates a column representing the gene type
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        return attributeDict[psl_lib.removeAlignmentNumber(aId)].geneType


class TranscriptType(object):
    """
    Creates a column representing the transcript type
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        return attributeDict[psl_lib.removeAlignmentNumber(aId)].transcriptType


class BeginStart(object):
    """

    Reports the 0-based position of the start codon in the transcript.

    Writes -1 to the database if this transcript either:
    1) has no thick region (thickStart == thickStop)
    2) has a thick region smaller than 3
    3) has no start codon in the thick window

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        s = t.getCds(seqDict)
        #ATG is the only start codon
        if len(s) == 0 or s[:3] != "ATG":
           return -1
        else:
            return t.cdsCoordinateToTranscript(0)


class CdsGap(object):
    """

    Are any of the CDS introns too short?

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict, 
                short_intron_size = 30):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) <= short_intron_size:
                    return 1
        return 0


class CdsMult3Gap(object):
    """

    Are any of the short CDS introns a multiple of 3?

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict, 
                short_intron_size = 30):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) <= short_intron_size and len(t.intronIntervals[i]) % 3 == 0:
                    return 1
        return 0


class CdsNonCanonSplice(object):
    """

    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def bad_splice(self, donor, acceptor):
        m = {"GT":"AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict, 
                short_intron_size = 30):
        if aId not in transcriptDict:
            return None              
        t = transcriptDict[aId]
        chromSeq = seqDict[t.chromosomeInterval.chromosome]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) > short_intron_size:
                    donor = chromSeq[t.intronIntervals[i].start : t.intronIntervals[i].start + 2]
                    acceptor = chromSeq[t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                    if self.bad_splice(donor, acceptor) is True:
                        return 1
        return 0


class CdsUnknownSplice(object):
    """

    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def bad_splice(self, donor, acceptor):
        m = {"GT":"AG", "GC":"AG", "AT":"AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict, 
                short_intron_size = 30):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        chromSeq = seqDict[t.chromosomeInterval.chromosome]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) > short_intron_size:
                    donor = chromSeq[t.intronIntervals[i].start : t.intronIntervals[i].start + 2]
                    acceptor = chromSeq[t.intronIntervals[i].stop - 2 : t.intronIntervals[i].stop]
                    if self.bad_splice(donor, acceptor) is True:
                        return 1
        return 0


class EndStop(object):
    """

    Looks at the end of the coding region (thickEnd) and sees if the last
    three bases are a stop codon ('TAA', 'TGA', 'TAG')

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None
        t = transcriptDict[aId]
        s = t.getProteinSequence(seqDict)
        if len(s) > 0 and s[-1] != "*":
            return 1
        return 0


class InFrameStop(object):
    """

    Looks for stop codons in frame in the coding sequence.

    Reports on in frame stop codons for each transcript.

    Records the 0-based transcript coordinate of an in-frame stop codon
    if it exists. Otherwise, records -1.

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None        
        t = transcriptDict[aId]
        #make sure this transcript has CDS
        #and more than 2 codons - can't have in frame stop without that
        cds_size = t.getCdsLength()
        if cds_size >= 9:
            for i in xrange(3, cds_size - 3, 3):
                c = t.cdsCoordinateToAminoAcid(i, seqDict)
                if c == "*":
                    return i
        return -1


class NoCds(object):
    """

    Looks to see if this transcript actually has a CDS

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None
        t = transcriptDict[aId]
        if t.getCdsLength() < 3:
            return 1
        return 0


class NumberScaffoldGap(object):
    """

    How many 100bp N runs are there in the alignment?
    100bp N runs are markers of scaffold gaps.

    """

    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        aln = alignmentDict[aId]
        destSeq = seqDict[aln.tName][aln.tStart : aln.tEnd].upper()
        if re.search("[N]{100}", destSeq) is not None:
            return 1
        return 0


class SourceChrom(object):
    """
    Creates a column representing the source chromosome
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        tName = psl_lib.removeAlignmentNumber(aId)
        if tName not in annotationDict:
            return None
        else:
            return annotationDict[tName].chromosomeInterval.chromosome


class SourceStart(object):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        tName = psl_lib.removeAlignmentNumber(aId)
        if tName not in annotationDict:
            return None
        else:
            return annotationDict[tName].chromosomeInterval.start


class SourceStop(object):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        tName = psl_lib.removeAlignmentNumber(aId)
        if tName not in annotationDict:
            return None
        else:
            return annotationDict[tName].chromosomeInterval.stop 


class SourceStrand(object):
    """
    Creates a column representing the source genomic strand.
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        tName = psl_lib.removeAlignmentNumber(aId)
        if tName not in annotationDict:
            s = None
        else:
            s = annotationDict[tName].chromosomeInterval.strand
        return seq_lib.convertStrand(s)


class DestChrom(object):
    """
    Creates a column representing the dest chromosome
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        if aId not in transcriptDict:
            return None
        else:
            return transcriptDict[aId].chromosomeInterval.chromosome


class DestStart(object):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        if aId not in transcriptDict:
            return None
        else:
            return transcriptDict[aId].chromosomeInterval.start


class DestStop(object):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        if aId not in transcriptDict:
            return None
        else:
            return transcriptDict[aId].chromosomeInterval.stop


class DestStrand(object):
    """
    Creates a column representing the dest genomic strand.
    """
    @staticmethod
    def _getType():
        return "TEXT"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):   
        if aId not in transcriptDict:
            s = None
        else:
            s = transcriptDict[aId].chromosomeInterval.strand
        return seq_lib.convertStrand(s)


class UnknownBases(object):
    """
    Counts the number of Ns in the target sequence within alignment blocks

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict):
        count = 0
        aln = alignmentDict[aId]
        if aln.strand == "+":
            for tStart, blockSize in izip(aln.tStarts, aln.blockSizes):
                seq = seqDict[aln.tName][tStart : tStart + blockSize]
                count += seq.count("N")
        else:
            #on negative strand the tStarts are (+) strand but blockSizes are in
            #transcript orientation
            for tStart, blockSize in izip(aln.tStarts, reversed(aln.blockSizes)):
                seq = seqDict[aln.tName][tStart : tStart + blockSize]
                count += seq.count("N")
        return count


class UtrGap(object):
    """

    Are any UTR introns too short?

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict,
                short_intron_size=30):
        if aId not in transcriptDict:
            return None      
        t = transcriptDict[aId]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is False and t.exons[i+1].containsCds() is False:
                if len(t.intronIntervals[i]) <= short_intron_size:
                    return 1
        return 0


class UtrNonCanonSplice(object):
    """

    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def bad_splice(self, donor, acceptor):
        m = {"GT":"AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict,
                short_intron_size=30):
        if aId not in transcriptDict:
            return None    
        t = transcriptDict[aId]
        chromSeq = seqDict[t.chromosomeInterval.chromosome]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is False and t.exons[i+1].containsCds() is False:
                if len(t.intronIntervals[i]) > short_intron_size:
                    donor = chromSeq[t.intronIntervals[i].start : t.intronIntervals[i].start + 2]
                    acceptor = chromSeq[t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                    if self.bad_splice(donor, acceptor) is True:
                        return 1
        return 0


class UtrUnknownSplice(object):
    """

    Are any of the UTR introns splice sites not of the
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def _getType():
        return "INTEGER"

    def bad_splice(self, donor, acceptor):
        m = {"GT":"AG", "GC":"AG", "AT":"AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def classify(self, aId, alignmentDict, refSeqDict, seqDict, attributeDict, transcriptDict, annotationDict,
                short_intron_size=30):
        if aId not in transcriptDict:
            return None        
        t = transcriptDict[aId]
        chromSeq = seqDict[t.chromosomeInterval.chromosome]
        for i in xrange(len(t.intronIntervals)):
            if t.exons[i].containsCds() is False and t.exons[i+1].containsCds() is False:
                if len(t.intronIntervals[i]) > short_intron_size:
                    donor = chromSeq[t.intronIntervals[i].start : t.intronIntervals[i].start + 2]
                    acceptor = chromSeq[t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                    if self.bad_splice(donor, acceptor) is True:
                        return 1
        return 0