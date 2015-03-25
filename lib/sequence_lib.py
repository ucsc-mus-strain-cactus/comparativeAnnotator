"""
Convenience library for sequence information, including 2bit and bed files.

fasta functionality was removed in favor of rapidly accessible 2bit.

Original Author: Dent Earl
Modified by Ian Fiddes
"""

import string
from itertools import izip
from collections import OrderedDict
from math import ceil, floor
from jobTree.src.bioio import fastaRead

#in case you are running the tests
try:
    from lib.twobit import TwoBitFile, TwoBitSequence
except:
    from twobit import TwoBitFile, TwoBitSequence

class Transcript(object):
    """
    Re-write this at some point, lol
    """
    
    __slots__ = ('chromosome', 'start', 'stop', 'name', 'score', 'strand', 'thickStart', 'thickStop', 'rgb', 
                 'blockCount', 'blockSizes', 'blockStarts', 'exonIntervals', 'intronIntervals', 'transcriptSize', 
                 'cdsSize', 'mRna', 'cds')
    
    def __init__(self, bed_tokens):
        self.chromosome = bed_tokens[0]
        self.start = int(bed_tokens[1])
        self.stop = int(bed_tokens[2])
        self.name = bed_tokens[3]
        self.score = int(bed_tokens[4])
        self.strand = convertStrand(bed_tokens[5])
        self.thickStart = int(bed_tokens[6])
        self.thickStop = int(bed_tokens[7])
        self.rgb = bed_tokens[8]
        self.blockCount = bed_tokens[9]
        self.blockSizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        self.blockStarts = [int(x) for x in bed_tokens[11].split(",") if x != ""]
        self.exonIntervals = self._getExonIntervals()
        self.intronIntervals = self._getIntronIntervals()

        self.transcriptSize = sum(x.stop - x.start for x in self.exonIntervals)
        self.cdsSize = self._getCdsLength()

    def __len__(self):
        return self.transcriptSize

    def __eq__(self, other):
        return (self.chromosomeInterval == other.chromosomeInterval and
                        self.name == other.name and
                        self.exonIntervals == other.exonIntervals and
                        self.score == other.score and
                        self.thickStart == other.thickStart and
                        self.thickStop == other.thickStop and
                        self.start == other.start and 
                        self.stop == other.end)

    def __cmp__(self, transcript):
        return cmp((self.chromosomeInterval, self.name),
                    (transcript.chromosomeInterval, transcript.name))

    def hashkey(self):
        """
        Return a string to use as dict key.
        """
        return '%s_%s_%d_%d' % (self.name, self.chromosomeInterval.chromosome, 
                self.chromosomeInterval.start, self.chromosomeInterval.stop)

    def getBed(self, rgb=None, name=None, start_offset=None, stop_offset=None):
        """
        Returns this transcript as a BED record with optional changes to rgb and name.
        If start_offset or stop_offset are set (chromosome coordinates), then this record will be changed to only 
        show results within that region, which is defined in chromosome coordinates.
        """
        if start_offset is not None and stop_offset is not None:
            assert start_offset <= stop_offset
        if start_offset is not None:
            assert start_offset >= self.start
        if stop_offset is not None:
            assert stop_offset <= self.stop
        if rgb is None:
            rgb = self.rgb
        if name is not None:
            name += "/" + self.name
        else:
            name = self.name
        if start_offset is None and stop_offset is None:
            return [self.chromosome, self.start, self.stop, name, self.score, convertStrand(self.strand), 
                    self.thickStart, self.thickStop, rgb, self.blockCount, self.blockSizes, self.blockStarts]
        elif start_offset == stop_offset:
            # no intron records
            assert self.chromosomeCoordinateToTranscript(start_offset) is not None 
            return [self.chromosome, start_offset, stop_offset, name, self.score, convertStrand(self.strand), 
                    start_offset, stop_offset, rgb, 1, 0, 0]
        
        def _moveStart(exonIntervals, blockCount, blockStarts, blockSizes, start, start_offset):
            toRemove = len([x for x in exonIntervals if x.start <= start_offset and x.stop <= start_offset])
            assert toRemove < len(exonIntervals)
            if toRemove > 0:
                blockCount -= toRemove
                blockSizes = blockSizes[toRemove:]
                start += blockStarts[toRemove]
                new_blockStarts = [0]
                for i in xrange(toRemove, len(blockStarts) - 1):
                    new_blockStarts.append(blockStarts[i + 1] - blockStarts[i] + new_blockStarts[-1])
                blockStarts = new_blockStarts
            if start_offset > start:
                blockSizes[0] += start - start_offset
                blockStarts[1:] = [x + start - start_offset for x in blockStarts[1:]]
                start = start_offset
            return start, blockCount, blockStarts, blockSizes
        
        def _moveStop(exonIntervals, blockCount, blockStarts, blockSizes, stop, start, stop_offset):
            toRemove = len([x for x in exonIntervals if x.stop >= stop_offset and x.start >= stop_offset])
            assert toRemove < len(exonIntervals)
            if toRemove > 0:
                blockCount -= toRemove
                blockSizes = blockSizes[:-toRemove]
                blockStarts = blockStarts[:-toRemove]
                assert len(blockSizes) == len(blockStarts)
                if len(blockSizes) == 0:
                    blockSizes = blockStarts = [0]
                    blockCount = 1
                stop = start + blockSizes[-1] + blockStarts[-1]
            if stop_offset < stop and stop_offset > start + blockStarts[-1]:
                blockSizes[-1] = stop_offset - start - blockStarts[-1] 
                stop = stop_offset
            return stop, blockCount, blockStarts, blockSizes
        
        blockCount = int(self.blockCount)
        blockStarts = self.blockStarts
        blockSizes = self.blockSizes
        start = self.start
        stop = self.stop
        thickStart = self.thickStart
        thickStop = self.thickStop
        
        if start_offset is not None and start_offset > start:
            start, blockCount, blockStarts, blockSizes = _moveStart(self.exonIntervals, blockCount, blockStarts, 
                                                                    blockSizes, start, start_offset)
        if stop_offset is not None and stop_offset < stop:
            stop, blockCount, blockStarts, blockSizes = _moveStop(self.exonIntervals, blockCount, blockStarts, 
                                                                  blockSizes, stop, start, stop_offset)
        if start > thickStart:
            thickStart = start
        if stop < thickStop:
            thickStop = stop
        if (start > thickStop and stop > thickStop) or (start < thickStart and stop < thickStart):
            thickStart = 0
            thickStop = 0
        blockStarts = ",".join(map(str, blockStarts))
        blockSizes = ",".join(map(str, blockSizes))
        return [self.chromosome, start, stop, name, self.score, convertStrand(self.strand), thickStart, thickStop, rgb,
                blockCount, blockSizes, blockStarts]

    def _getExonIntervals(self):
        """
        Gets a list of exon intervals in chromosome coordinate space.
        These exons are on (+) strand ordering regardless of transcript strand.
        This means (-) strand genes will be represented backwards
        """
        exons = []
        for block_size, block_start in izip(self.blockSizes, self.blockStarts):
            exons.append(ChromosomeInterval(self.chromosome, self.start + block_start, self.start + block_start + \
                         block_size, self.strand))
        return exons

    def _getIntronIntervals(self):
        """
        Get a list of ChromosomeIntervals representing the introns for this
        transcript. The introns are in *+ strand of CHROMOSOME* ordering,
        not the order that they appear in the transcript!
        """
        introns = []
        prevExon = None
        for exon in self.exonIntervals:
            if prevExon is not None:
                assert exon.start > prevExon.stop
                assert exon.strand == prevExon.strand
                introns.append(ChromosomeInterval(self.chromosome, prevExon.stop, exon.start, exon.strand))
            prevExon = exon
        return introns

    def _getCdsLength(self):
        """
        Returns the length of the CDS.
        """
        l = 0
        for e in self.exonIntervals:
            if self.thickStart < e.start and e.stop < self.thickStop:
                # squarely in the CDS
                l += e.stop - e.start
            elif (e.start <= self.thickStart and e.stop < self.thickStop
                        and self.thickStart < e.stop):
                # thickStart marks the start of the CDS
                l += e.stop - self.thickStart
            elif e.start <= self.thickStart and self.thickStop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                l += self.thickStop - self.thickStart
            elif (self.thickStart < e.start and self.thickStop <= e.stop
                        and e.start < self.thickStop):
                # thickStop marks the end of the CDS
                l += self.thickStop - e.start
        return l

    def getCdsLength(self):
        return self.cdsSize

    def getTranscriptLength(self):
        return self.transcriptSize

    def hasCds(self):
        if self.cdsSize > 0:
            return True
        else:
            return False

    def getMRna(self, seqDict):
        """
        Returns the mRNA sequence for this transcript based on a TwoBitFile object.
        and the start/end positions and the exons. Sequence returned in
        5'-3' transcript orientation.
        """
        if hasattr(self, "mRna"):
            return self.mRna
        sequence = seqDict[self.chromosome]
        assert self.stop <= len(sequence)
        s = []
        for e in self.exonIntervals:
            s.append(sequence[e.start : e.stop])
        if self.strand is True:
            mRna = "".join(s)
        else:
            mRna = reverseComplement("".join(s))
        self.mRna = mRna.upper()
        return mRna

    def getSequence(self, seqDict):
        """
        Returns the entire chromosome sequence for this transcript, (+) strand orientation.
        """
        sequence = seqDict[self.chromosomeInterval.chromosome]
        return sequence[self.start:self.stop]

    def getCds(self, seqDict):
        """
        Return the CDS sequence (as a string) for the transcript 
        (based on the exons) using a sequenceDict as the sequence source.
        The returned sequence is in the correct 5'-3' orientation (i.e. it has
        been reverse complemented if necessary).
        """
        if hasattr(self, "cds"):
            return self.cds
        sequence = seqDict[self.chromosome]
        assert self.stop <= len(sequence)
        #make sure this isn't a non-coding gene
        if self.thickStart == self.thickStop == 0:
            return ""
        s = []
        for e in self.exonIntervals:
            if self.thickStart < e.start and e.stop < self.thickStop:
                # squarely in the CDS
                s.append(sequence[e.start : e.stop])
            elif (e.start <= self.thickStart and e.stop < self.thickStop
                        and self.thickStart < e.stop):
                # thickStart marks the start of the CDS
                s.append(sequence[self.thickStart : e.stop])
            elif e.start <= self.thickStart and self.thickStop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                s.append(sequence[self.thickStart : self.thickStop])
            elif (self.thickStart < e.start and self.thickStop <= e.stop
                        and e.start < self.thickStop):
                # thickStop marks the end of the CDS
                s.append(sequence[e.start : self.thickStop])
        if not self.strand:
            cds = reverseComplement("".join(s)).upper()
        else:
            cds = "".join(s).upper()
        self.cds = cds
        return cds

    def getProteinSequence(self, seqDict):
        """
        Returns the translated protein sequence for this transcript in single
        character space.
        """
        cds = self.getCds(seqDict)
        if len(cds) < 3:
            return ""
        return translateSequence(self.getCds(seqDict))

    def intronSequenceIterator(self, seqDict):
        """
        Iterates over intron sequences in transcript order and strand
        """
        if self.strand is True:
            for intron in self.intronIntervals:
                yield intron.getSequence(seqDict, strand=True)
        else:
            for intron in reversed(self.intronIntervals):
                yield intron.getSequence(seqDict, strand=True)

    def getIntronSequences(self, seqDict):
        """
        Wrapper for intronSequenceIterator that returns a list of sequences.
        """
        return list(self.intronSequenceIterator(seqDict))

    def transcriptCoordinateToCds(self, p):
        if p is None: return None
        # find the thickStart, thickStop offsets in exon coordinates
        exonThickStart, exonthickStop = None, None
        x = 0  # exon coordinate
        for e in self.exonIntervals:
          length = e.stop - e.start
          if exonThickStart is None and e.start >= self.thickStart:
            # thickStart fell between exons
            exonThickStart = x
          if exonThickStart is None and e.stop > self.thickStart:
            # exon contains thickStart
            exonThickStart = x + self.thickStart - e.start
          if exonthickStop is None and e.start >= self.thickStop:
            # thickStop fell between exons
            exonthickStop = x
          if exonthickStop is None and e.stop >= self.thickStop:
            # exon contains thickStop
            exonthickStop = x + self.thickStop - e.start
          x += length
        if self.strand is False:
          exonThickStart, exonthickStop = exonthickStop, exonThickStart
          exonThickStart = x - exonThickStart
          exonthickStop = x - exonthickStop
        if p < exonThickStart:
          return None
        if p >= exonthickStop:
          return None
        return p - exonThickStart


    def transcriptCoordinateToChromosome(self, p):
        if p is None: return None
        if p < 0:
          return None
        if p >= self.transcriptSize:
          return None
        c = 0  # cumulative position through exon space
        if self.strand is not True:
          p = self.transcriptSize - 1 - p
        e_start = self.exonIntervals[0].start
        for e in self.exonIntervals:
          if p < c + e.stop - e.start:
            # the position is within this exon
            return p - c + e.start
          else:
            # sorry mario, your position is in another exon
            c += e.stop - e.start
        assert False   # we should never get here

    def cdsCoordinateToTranscript(self, p):
        if p is None: return None
        if p < 0: return None
        if p >= self.cdsSize: return None
        if self.strand is True:
          # positive strand, offset is first exon start to thickStart
          for e in self.exonIntervals:
            if e.start < self.thickStart and e.stop <= self.thickStart:
              # add the whole exon to the offset
              p += e.stop - e.start
            elif e.start < self.thickStart and self.thickStart <= e.stop:
              # only add the thin part of this exon
              p += self.thickStart - e.start
              break
        else:
          for e in reversed(self.exonIntervals):
            if self.thickStop < e.start and self.thickStop < e.stop:
              # add the whole exon to the offset
              p += e.stop - e.start
            elif e.start < self.thickStop and self.thickStop < e.stop:
              # only add the thin part of this exon
              p += e.stop -  self.thickStop
              break
        return p

    def cdsCoordinateToChromosome(self, p):
        if p is None: return None
        if p < 0: return None
        if p >= self.cdsSize: return None
        q = self.cdsCoordinateToTranscript(p)
        if q >= self.transcriptSize: return None
        return self.transcriptCoordinateToChromosome(q)

    def chromosomeCoordinateToTranscript(self, p):
        if p is None: return None
        if self.strand is True:
          def _stranded(v): return v
        else:
          def _stranded(v):
            return self.transcriptSize - 1 - v
        c = 0  # cumulative position through exon space
        e_start = self.exonIntervals[0].start
        for e in self.exonIntervals:
          if p < e.start:
            # p is not in an exon
            return None
          if p < e.stop:
            # the position is within this exon
            return _stranded(c + p - e.start)
          else:
            # sorry mario, your position is in another exon
            c += e.stop - e.start
        return None

    def chromosomeCoordinateToCds(self, p):
        if p is None: return None
        if p < 0:
          return None
        if p >= self.stop:
          return None
        q = self.chromosomeCoordinateToTranscript(p)
        if q is None:
          return None
        if q < 0:
          return None
        if q >= self.transcriptSize:
          return None
        return self.transcriptCoordinateToCds(q)

    def cdsCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a CDS-relative position and a sequenceDict object that contains this
        transcript and returns the amino acid at that CDS position.
        Returns None if this is invalid.
        """
        cds = self.getCds(seqDict)
        if p is None or p >= len(cds) or p < 0:
            return None
        #we add 0.1 to the ceiling to make multiples of 3 work
        start, stop = int(floor(p / 3.0) * 3),  int(ceil((p + 0.1) / 3.0) * 3)
        if stop - start != 3:
            return None
        codon = cds[start:stop]
        return codonToAminoAcid(codon)

    def transcriptCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a transcript coordinate position and a sequenceDict object that contains
        this transcript and returns the amino acid at that transcript position.
        If this position is not inside the CDS returns None.
        """
        cds_pos = self.transcriptCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, seqDict)

    def chromosomeCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a chromosome coordinate and a sequenceDict object that contains this
        transcript and returns the amino acid at that chromosome position.
        Returns None if this position is not inside the CDS.
        """
        cds_pos = self.chromosomeCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, seqDict)


class ChromosomeInterval(object):
    """
    Represents an interval of a chromosome. BED coordinates, strand is True,
    False or None (if no strand)
    """
    
    __slots__ = ('chromosome', 'start', 'stop', 'strand')    # conserve memory
    
    def __init__(self, chromosome, start, stop, strand):
        self.chromosome = str(chromosome)
        self.start = int(start)    # 0 based
        self.stop = int(stop)    # exclusive
        assert(strand in [True, False, None])
        self.strand = strand    # True or False

    def __eq__(self, other):
        return (self.chromosome == other.chromosome and self.start == other.start and  self.stop == other.stop and
                self.strand == other.strand)

    def __cmp__(self, cI):
        return cmp((self.chromosome, self.start, self.stop, self.strand), (cI.chromosome, cI.start, cI.stop, cI.strand))

    def __len__(self):
        return self.size()

    def contains(self, other):
        """ Check the other chromosomeInterval to see if it is contained by this
        CI. If it is not contained return False, else return True.
        """
        if not isinstance(other, ChromosomeInterval):
            raise RuntimeError('ChromosomeInterval:contains expects '
                                'ChromosomeInterval, not %s' % other.__class__)
        if self.chromosome != other.chromosome:
            return False
            # self  |----*
            # other         *----|
        if self.stop <= other.start:
            return False
            # self          *----|
            # other |----*
        if self.start >= other.stop:
            return False
            # self    *------|
            # other *----|
        if self.start > other.start:
            return False
            # self  |-----*
            # other    |----*
        if self.stop < other.stop:
            return False
        return True

    def size(self):
        return self.stop - self.start

    def getBed(self, name):
        """
        Returns BED tokens representing this interval. Requires a name.
        """
        return [self.chromosome, self.start, self.stop, name, 0, convertStrand(self.strand)]

    def getSequence(self, seqDict, strand=True):
        """
        Returns a string representing the sequence of this interval. If strand is True, converts strand.
        """
        seq = seqDict[self.chromosome][self.start:self.stop]
        if self.strand is True:
            return seq
        else:
            return reverseComplement(seq)


class Attribute(object):
    """
    Stores attributes from the gencode attribute file.
    """
    
    __slots__ = ("geneId", "geneName", "geneType", "transcriptId", "transcriptType")
    
    def __init__(self, geneId, geneName, geneType, transcriptId, transcriptType):
        self.geneId = geneId
        self.geneName = geneName
        self.geneType = geneType
        self.transcriptId = transcriptId
        self.transcriptType = transcriptType


def convertStrand(s):
    """
    Given a potential strand value, converts either from True/False/None
    to +/-/None depending on value
    """
    assert s in [True, False, None, "+", "-"]
    if s == True: return "+"
    elif s == False: return "-"
    elif s == None: return None
    elif s == "-": return False
    elif s == "+": return True


_complement = string.maketrans("ATGC","TACG")

def complement(seq):
  """
  given a sequence, return the complement.
  """
  return seq.translate(_complement)


def reverseComplement(seq):
  """
  Given a sequence, return the reverse complement.
  """
  return seq.translate(_complement)[::-1]


_codonTable = {
    'ATG': 'M',
    'TAA': '*', 'TAG': '*', 'TGA': '*', 'TAR': '*', 'TRA': '*',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCN': 'A',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R','AGA': 'R',
    'AGG': 'R', 'CGN': 'R', 'MGR': 'R',
    'AAT': 'N', 'AAC': 'N', 'AAY': 'N',
    'GAT': 'N', 'GAC': 'N', 'GAY': 'N',
    'TGT': 'C', 'TGC': 'C', 'TGY': 'C',
    'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q',
    'GAA': 'E', 'GAG': 'E', 'GAR': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
    'CAT': 'H', 'CAC': 'H', 'CAY': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATH': 'I',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
    'CTG': 'L', 'YTR': 'L', 'CTN': 'L',
    'AAA': 'K', 'AAG': 'K', 'AAR': 'K',
    'TTT': 'F', 'TTC': 'F', 'TTY': 'F',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S',
    'AGC': 'S', 'TCN': 'S', 'AGY': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACN': 'T',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTN': 'V',
    }

def codonToAminoAcid(c):
    """
    Given a codon C, return an amino acid or ??? if codon unrecognized.
    Codons could be unrecognized due to ambiguity in IUPAC characters.
    """
    if c is None: return None
    c = c.upper()
    if c in _codonTable:
        return _codonTable[c]
    return '?'


def translateSequence(sequence):
    """
    Translates a given DNA sequence to single-letter amino acid
    space. If the sequence is not a multiple of 3 it will be truncated
    silently.
    """
    result = []
    for i in xrange(0, len(sequence)- len(sequence) % 3, 3):
        result.append(codonToAminoAcid(sequence[i : i + 3]))
    return "".join(result)


def readCodons(seq):
    """
    Provides an iterator that reads through a sequence one codon at a time.
    """
    l = len(seq)
    for i in xrange(0, l, 3):
        if i + 3 <= l:
            yield seq[i:i + 3]

def readCodonsWithPosition(seq):
    """
    Provides an iterator that reads through a sequence one codon at a time,
    returning both the codon and the start position in the sequence.
    """
    l = len(seq)
    for i in xrange(0, l, 3):
        if i + 3 <= l:
            yield i, seq[i:i + 3]

def readTwoBit(file_path):
    """
    Returns a dictionary that can randomly access two bit files.
    Acts as a wrapper around the TwoBitFile class in twobitreader.py.
    """
    return TwoBitFile(file_path)


def getSequenceDict(file_path):
    """
    Returns a dictionary of fasta records.
    """
    fastaDict = {}
    for name, seq in fastaRead(file_path):
        fastaDict[name] = seq
    return fastaDict


def getTranscripts(bedFile):
    """
    Given a path to a standard BED file and a details BED, return a list of
    Transcript objects.
    """
    transcripts = []
    bedFile = open(bedFile, 'r')
    for t in transcriptIterator(bedFile):
        transcripts.append(t)
    return transcripts


def transcriptListToDict(transcripts, noDuplicates=False):
    """
    Given a list af Transcript objects, attempt to transform them into a dict
    of lists. key is transcript name, value is list of Transcript objects.
    If NODUPLICATES is true, then the value will be a single Transcript object.
    """
    result = {}
    for t in transcripts:
        if t.name not in result:
            result[t.name] = []
        else:
            if noDuplicates:
                raise RuntimeError('transcriptListToDict: Discovered a duplicate transcript %s %s' 
                                   % (t.name, t.chromosome))
        if noDuplicates:
            result[t.name] = t
        else:
            result[t.name].append(t)
    return result


def tokenizeBedStream(bedStream):
    """
    Iterator through bed file, returning lines as list of tokens
    """
    for line in bedStream:
        if line != '':
            tokens = line.split()
            yield tokens


def tokenizeGenePredStream(genePredStream):
    """
    Iterator through gene pred file, returning lines as list of tokens
    """
    for line in bedStream:
        if line != '':
            tokens = line.split("\t")
            tokens[-1].rstrip()
            yield tokens


def transcriptIterator(transcriptsBedStream):
    """
    Iterates over the transcript detailed in the bed stream producing Transcript objects.
    """
    for tokens in tokenizeBedStream(transcriptsBedStream):
        yield Transcript(tokens)


def getTranscriptAttributeDict(attributeFile):
    """
    Returns a dictionary mapping the transcript ID to an Attribute object.
    This stores all of the relevant information from the gencode attributes file.
    """
    attribute_dict = {}
    with open(attributeFile) as f: 
        for line in f:
            line = line.split("\t")
            if line[0] == "geneId": 
                continue
            geneID, geneName, geneType, geneStatus, transcriptID, transcriptName, transcriptType, transcriptStatus, \
                    havanaGeneID, havanaTranscriptID, ccdsID, level, transcriptClass = line
            attribute_dict[transcriptID] = Attribute(geneID, geneName, geneType, transcriptID, transcriptType)
    return attribute_dict


def intervalToBed(t, interval, rgb, name):
    """
    If you are turning interval objects into BED records, look here. t is a transcript object.
    Interval objects should always have start <= stop (+ strand chromosome ordering)
    """
    assert interval.stop >= interval.start
    return [interval.chromosome, interval.start, interval.stop, name + "/" + t.name, 0, convertStrand(interval.strand), 
            interval.start, interval.stop, rgb, 1, interval.stop - interval.start, 0]


def spliceIntronIntervalToBed(t, intronInterval, rgb, name):
    """
    Specific case of turning an intron interval into the first and last two bases (splice sites)
    """
    interval = intronInterval
    assert interval.stop >= interval.start
    blockStarts = "0,{}".format(interval.stop - 2)
    return [inteval.chromosome, interval.start, interval.stop, name + "/" + t.name, 0, convertStrand(interval.strand), 
            interval.start, interval.stop, rgb, 2, "2,2", blockStarts]


def transcriptToBed(t, rgb, name):
    """
    Convenience function for pulling BED tokens from a Transcript object.
    """
    return t.getBed(rgb, name)


def transcriptCoordinateToBed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in TRANSCRIPT coordinate space and returns
    a list in BED format with the specified RGB string (128,0,0 or etc) and name.
    """
    chromStart = t.transcriptCoordinateToChromosome(start)
    chromStop = t.transcriptCoordinateToChromosome(stop)
    return chromosomeCoordinateToBed(t, chromStart, chromStop, rgb, name)


def cdsCoordinateToBed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CDS coordinate space and returns
    a list in BED format with the specified RGB string (128,0,0 or etc) and name.
    """
    if t.strand is False:
        chromStart = t.cdsCoordinateToChromosome(stop)
        chromStop = t.cdsCoordinateToChromosome(start) + 1
    else:
        chromStart = t.cdsCoordinateToChromosome(start)
        chromStop = t.cdsCoordinateToChromosome(stop) + 1
    return chromosomeCoordinateToBed(t, chromStart, chromStop, rgb, name)


def chromosomeCoordinateToBed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CHROMOSOME coordinate space and returns
    a list in BED format with the specified RGB string and name.
    """
    assert start != None and stop != None, (t.name, start, stop, name)
    assert stop >= start, (t.name, start, stop, name)
    return t.getBed(name=name, rgb=rgb, start_offset=start, stop_offset=stop)


def chromosomeRegionToBed(t, start, stop, rgb, name):
    """
    This is different from chromosomeCoordinateToBed - this function will not resize the BED information
    for the input transcript, but instead be any coordinate on the chromosome.
    """
    assert start != None and stop != None, (t.name, start, stop, name)
    assert stop >= start, (t.name, start, stop, name)
    return [t.chromosome, start, stop, name + "/" + t.name, 0, t.strand, start, stop, rgb, 1, stop - start, 0]        
