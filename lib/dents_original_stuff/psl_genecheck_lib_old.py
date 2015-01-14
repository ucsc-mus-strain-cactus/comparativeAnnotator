"""
THIS IS THE ORIGINAL LIBRARY FROM DENT

convenience library for parsing genecheck psl files
author: dent earl
"""
from argparse import ArgumentTypeError
from collections import deque
import os
import sys


class Sequence(object):
  """ Represents a sequence of DNA.
  """
  __slots__ = ('name', '_sequence', '_length')  # conserve memory
  def __init__(self, name, sequence):
    self.name = name  # chromosome or scaffold name
    self._sequence = sequence  # ACGTs
    self._length = len(sequence)
  def setSequence(self, seq):
    self._sequence = seq
    self._length = len(seq)
  def getSequence(self):
    return self._sequence
  def getLength(self):
    return self._length
  def setUpper(self):
    self._sequence = self._sequence.upper()
  def getNucleotide(self, pos, relativeStrand=True, complementNuc=False):
    """ return the single nucleotide that resides at 0-based position POS.
    """
    if relativeStrand:
      n = self.sliceSequence(pos, pos + 1)
      if complementNuc:
        n = complement(n)
    else:
      # using slice sequence with '-' automatically complements
      n = self.sliceSequence(pos, pos + 1, relativeStrand='-')
    return n
  def sliceSequence(self, start, stop, relativeStrand='+'):
    """ return the proper slice of the sequence.
    BED format coordinates: 0 based start, stop is exclusive
    [start, stop). E.g. Sequence.sliceSequence(0, 3) returns a string length 3.
    """
    assert(start < stop)
    if relativeStrand == '+':
      return self._sequence[start:stop]
    elif relativeStrand == '-':
      # 0 1 2 3 4 5 6 7 8 9  +
      # 9 8 7 6 5 4 3 2 1 0  -
      #                  + strand | - strand
      #               |    [7, 8) = [2, 3)
      #   |-----|          [1, 5) = [5, 9)
      # |---------|        [0, 6) = [4, 10)
      #       |---------|  [3, 9) = [1, 7)
      a = self._length - stop
      b = self._length - stop + (stop - start)
      return reverseComplement(self._sequence[a:b])
    else:
      raise RuntimeError('Unanticipated relativeStrand: %s'
                         % str(relativeStrand))


class PslRow(object):
  """ Represents a single row in a PSL file.
  http://genome.ucsc.edu/FAQ/FAQformat.html#format2
  """
  __slots__ = ('matches', 'misMatches', 'repMatches', 'nCount',
               'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert',
               'strand', 'qName', 'qSize', 'qStart', 'qEnd',
               'tName', 'tSize', 'tStart', 'tEnd', 'blockCount',
               'blockSizes', 'qStarts', 'tStarts')  # conserve memory
  def __init__(self, line):
    data = line.split()
    assert(len(data) == 21)
    self.matches = int(data[0])
    self.misMatches = int(data[1])
    self.repMatches = int(data[2])
    self.nCount = int(data[3])
    self.qNumInsert = int(data[4])
    self.qBaseInsert = int(data[5])
    self.tNumInsert = int(data[6])
    self.tBaseInsert = int(data[7])
    self.strand = data[8]
    self.qName = data[9]
    self.qSize = int(data[10])
    self.qStart = int(data[11])
    self.qEnd = int(data[12])
    self.tName = data[13]
    self.tSize = int(data[14])
    self.tStart = int(data[15])
    self.tEnd = int(data[16])
    self.blockCount = int(data[17])
    # lists of ints
    self.blockSizes = [int(x) for x in data[18].split(',') if x]
    self.qStarts = [int(x) for x in data[19].split(',') if x]
    self.tStarts = [int(x) for x in data[20].split(',') if x]
  def hashkey(self):
    """ return a string to use as dict key.
    """
    return '%s_%s_%d_%d' % (self.qName, self.tName, self.tStart, self.tEnd)
  def targetCoordinateToQuery(self, p):
    """ Take position P in target coordinates (positive) and convert it
    to query coordinates (positive). If P is not in target coordinates throw
    assert, if P does not map to query coordinates return None.
    """
    if p < self.tStart: return None
    if p >= self.tEnd: return None
    if self.strand not in ['+', '-']:
      raise RuntimeError('Unanticipated strand: %s' % self.strand)
    for i, t in enumerate(self.tStarts):
      if p < t:
        continue
      if p >= t + self.blockSizes[i]:
        continue
      # p must be in block
      offset = p - t
      if self.strand == '+':
        return self.qStarts[i] + offset
      else:
        return self.qSize - (self.qStarts[i] + offset) - 1
    return None
  def queryCoordinateToTarget(self, p):
    """ Take position P in query coordinates (positive) and convert it
    to target coordinates (positive). If P is not in query coordinates throw
    assert, if P does not map to target coordinates return None.
    """
    # this is the easier one to write
    if self.strand == '+':
      pass
    elif self.strand == '-':
      p = self.qSize - p - 1
    else:
      raise RuntimeError('Unanticipated strand: %s' % self.strand)
    if p < self.qStart: return None
    if p >= self.qEnd: return None
    for i, q in enumerate(self.qStarts):
      if p < q:
        continue
      if p >= q + self.blockSizes[i]:
        continue
      # p must be in block
      offset = p - q
      return self.tStarts[i] + offset
    return None
  def pslString(self):
    """ return SELF as a psl formatted line.
    """
    s = ('%d %d %d %d %d %d %d %d %s %s %d %d %d %s %d %d %d %d %s %s %s' %
         (self.matches,
          self.misMatches,
          self.repMatches,
          self.nCount,
          self.qNumInsert,
          self.qBaseInsert,
          self.tNumInsert,
          self.tBaseInsert,
          self.strand,
          self.qName,
          self.qSize,
          self.qStart,
          self.qEnd,
          self.tName,
          self.tSize,
          self.tStart,
          self.tEnd,
          self.blockCount,
          # lists of ints
          ','.join([str(b) for b in self.blockSizes]),
          ','.join([str(b) for b in self.qStarts]),
          ','.join([str(b) for b in self.tStarts])))
    return s


"""The following data types are used for iterating over gene-check-detail and
gene-check bed files.
An example of entries from such files:

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check-details.bed
1       2812370 2812372 noStop/ENSMUST00000065527.4
1       2812370 2812372 noStop/ENSMUST00000095795.4

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check.bed
1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729
128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301
050,301352
1       2812346 3113783 ENSMUST00000095795.4    0       -       2812370 3038729
128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301
050,301352
"""

class ChromosomeInterval(object):
  """ Represents an interval of a chromosome. BED coordinates, strand is True,
  False or None (if no strand)
  """
  __slots__ = ('chromosome', 'start', 'stop', 'strand')  # conserve memory
  def __init__(self, chromosome, start, stop, strand):
    self.chromosome = str(chromosome)
    self.start = int(start)  # 0 based
    self.stop = int(stop)  # exclusive
    assert(strand in [True, False, None])
    self.strand = strand  # True or False

  def __eq__(self, other):
    return (self.chromosome == other.chromosome and
            self.start == other.start and
            self.stop == other.stop and
            self.strand == other.strand)

  def __cmp__(self, cI):
    return cmp((self.chromosome, self.start, self.stop, self.strand),
               (cI.chromosome, cI.start, cI.stop, cI.strand))

  def contains(self, other):
    """ Check the other chromosomeInterval to see if it is contained by this
    CI. If it is not contained return False, else return True.
    """
    if not isinstance(other, ChromosomeInterval):
      raise RuntimeError('ChromosomeInterval:contains expects '
                         'ChromosomeInterval, not %s' % other.__class__)
    # print 'testing contains! me:[%d,%d), them:[%d,%d)' % (self.start, self.stop, other.start, other.stop)
    if self.chromosome != other.chromosome:
      # print 'nope'
      return False
    # self  |----*
    # other         *----|
    if self.stop <= other.start:
      # print 'nope'
      return False
    # self          *----|
    # other |----*
    if self.start >= other.stop:
      # print 'nope'
      return False
    # self    *------|
    # other *----|
    if self.start > other.start:
      # print 'nope'
      return False
    # self  |-----*
    # other    |----*
    if self.stop < other.stop:
      #print 'nope'
      return False
    # print 'yup!'
    return True

  def size(self):
    return self.stop - self.start


class TranscriptAnnotation(object):
  """ Represents an annotation of a transcript, from one of the
  classification bed files
  """
  # conserve memory with __slots__
  __slots__ = ('chromosomeInterval', 'name', 'labels', '_itemRgb')
  def __init__(self, chromosomeInterval, name, label):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.labels = deque(label)  # need to be able to append left, -> deque
    self._itemRgb = ''

  def addLabel(self, label, prepend=False):
    """ maintain the ordering of self.labels but prevent duplicates
    """
    if label not in self.labels:
      if prepend:
        self.labels.appendleft(label)
      else:
        self.labels.append(label)

  def bedString(self):
    strandChar = '-'
    if self.chromosomeInterval.strand:
      strandChar = '+'
    if self._itemRgb == '':
      return '\t'.join([self.chromosomeInterval.chromosome,
                        str(self.chromosomeInterval.start),
                        str(self.chromosomeInterval.stop),
                        '/'.join(list(self.labels) + [self.name])])
    else:
      return '\t'.join([self.chromosomeInterval.chromosome,
                        str(self.chromosomeInterval.start),
                        str(self.chromosomeInterval.stop),
                        '/'.join(list(self.labels) + [self.name]),
                        '0', strandChar,
                        str(self.chromosomeInterval.start),
                        str(self.chromosomeInterval.stop),
                        self._itemRgb])

  def __eq__(self, other):
    return (self.chromosomeInterval == other.chromosomeInterval and
            self.name == other.name and
            self.labels == other.labels)

  def __cmp__(self, annotation):
    """ Sort by chromosome interval, then name
    """
    return cmp((self.chromosomeInterval, self.name),
               (annotation.chromosomeInterval, annotation.name))

class Transcript(object):
  """ Represent a transcript and its annotations
  """
  __slots__ = ('chromosomeInterval', 'name', 'exons', 'annotations',
               'score', 'thickStart', 'thickEnd', 'itemRgb')  # conserve memory
  def __init__(self, chromosomeInterval, name, exons, annotations,
               score, thickStart, thickEnd, itemRgb):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.exons = exons  # list of chromosome intervals
    self.annotations = annotations  # list of TranscriptAnnotation() objects
    # Bed fields
    self.score = score
    self.thickStart = thickStart  # int
    self.thickEnd = thickEnd  # int
    self.itemRgb = itemRgb

  def __eq__(self, other):
    return (self.chromosomeInterval == other.chromosomeInterval and
            self.name == other.name and
            self.exons == other.exons and
            self.annotations == other.annotations and
            self.score == other.score and
            self.thickStart == other.thickStart and
            self.thickEnd == other.thickEnd and
            self.itemRgb == other.itemRgb)

  def hashkey(self, deuniquify=False):
    """ return a string to use as dict key.
    """
    if deuniquify:
      return '%s_%s_%d_%d' % (removeAlignmentNumber(self.name),
                              self.chromosomeInterval.chromosome,
                              self.chromosomeInterval.start,
                              self.chromosomeInterval.stop)
    else:
      return '%s_%s_%d_%d' % (self.name, self.chromosomeInterval.chromosome,
                              self.chromosomeInterval.start,
                              self.chromosomeInterval.stop)

  def mRnaCoordinateToCodon(self, p):
    """ Take position P with 0-based mRNA-relative position and convert it
    to 0-based (codon, codon position) tuple.
    """
    # (strand irrelevant)
    #       0    5    10
    #       |    |    |
    # mRNA  +++++++++++
    # codon +++++++++++
    #       |  |  |  |
    #       0,0   2,0
    if p is None: return None
    assert(p >= 0)
    assert(p < sum([(e.stop - e.start) for e in self.exons]))  # could be a tighter bound
    return (int(p / 3), p % 3)

  def codonCoordinateToMRna(self, p):
    """ Take (codon, codon position) tuple P and convert it
    to 0-based mRNA-relative position.
    """
    # (strand irrelevant)
    #       0    5    10
    #       |    |    |
    # mRNA  +++++++++++
    # codon +++++++++++
    #       |  |  |  |
    #       0,0   2,0
    if p is None: return None
    assert(p[0] >= 0)
    assert(p[1] >= 0)
    return p[0] * 3 + p[1]

  def codonCoordinateToChromosome(self, p):
    """ Take (codon, codon position) tuple P and convert it
    to 0-based chromosome-relative position.
    """
    m = self.codonCoordinateToMRna(p)
    return self.mRnaCoordinateToChromosome(m)

  def chromosomeCoordinateToCodon(self, p):
    """ Take 0-based chromosome-relative P and convert it
    to (codon, codon position) tuple position.
    """
    m = self.chromosomeCoordinateToMRna(p)
    return self.mRnaCoordinateToCodon(m)

  def mRnaCoordinateToExon(self, p):
    """ Take position P with 0-based mRNA-relative position and convert it
    to 0-based exon-relative position.
    """
    assert(len(self.exons))
    if p is None: return None
    if p < 0: return None
    if p >= sum([(e.stop - e.start) for e in self.exons]): return None
    # positive strand
    # chromosome +++++++++++
    #            |    |    |
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    # mrna          ++ ++++
    #               |     |
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    if self.chromosomeInterval.strand:
      # positive strand, offset is first exon start to thickStart
      for e in self.exons:
        if e.start < self.thickStart and e.stop <= self.thickStart:
          # add the whole exon to the offset
          p += e.stop - e.start
        elif e.start < self.thickStart and self.thickStart <= e.stop:
          # only add the thin part of this exon
          p += self.thickStart - e.start
          break
    else:
      for e in reversed(self.exons):
        if self.thickEnd < e.start and self.thickEnd < e.stop:
          # add the whole exon to the offset
          p += e.stop - e.start
        elif e.start < self.thickEnd and self.thickEnd < e.stop:
          # only add the thin part of this exon
          p += e.stop -  self.thickEnd
          break
    return p

  def exonCoordinateToMRna(self, p):
    """ Take position P with 0-based exon-relative position and convert it
    to 0-based mRNA-relative position. If position does not exist in mRNA,
    return None.
    """
    if p is None: return None
    # find the thickStart, thickEnd offsets in exon coordinates
    exonThickStart, exonThickEnd = None, None
    x = 0  # exon coordinate
    for e in self.exons:
      length = e.stop - e.start
      if exonThickStart is None and e.start >= self.thickStart:
        # thickStart fell between exons
        exonThickStart = x
      if exonThickStart is None and e.stop > self.thickStart:
        # exon contains thickStart
        exonThickStart = x + self.thickStart - e.start
      if exonThickEnd is None and e.start >= self.thickEnd:
        # thickEnd fell between exons
        exonThickEnd = x
      if exonThickEnd is None and e.stop >= self.thickEnd:
        # exon contains thickEnd
        exonThickEnd = x + self.thickEnd - e.start
      x += length
    if not self.chromosomeInterval.strand:
      exonThickStart, exonThickEnd = exonThickEnd, exonThickStart
      exonThickStart = x - exonThickStart
      exonThickEnd = x - exonThickEnd
    if p < exonThickStart:
      return None
    if p >= exonThickEnd:
      return None
    return p - exonThickStart

  def mRnaCoordinateToChromosome(self, p):
    """ Take position P with 0-based mRNA-relative position and convert it
    to 0-based chromosome-relative position.
    """
    assert(len(self.exons))
    if p is None: return None
    if p < 0: return None
    limit = sum([(e.stop - e.start) for e in self.exons])
    if p >= limit: return None
    p = self.mRnaCoordinateToExon(p)
    if p >= limit: return None
    return self.exonCoordinateToChromosome(p)

  def exonCoordinateToChromosome(self, p):
    """ Take position P with 0-based exon-relative position and convert it
    to 0-based chromosome-relative position.
    """
    if p is None: return None
    if p < 0:
      return None
    if p >= sum([(e.stop - e.start) for e in self.exons]):
      return None
    assert(len(self.exons))
    c = 0  # cumulative position through exon space
    if not self.chromosomeInterval.strand:
      p = sum([(e.stop - e.start) for e in self.exons]) - 1 - p
    e_start = self.exons[0].start
    for e in self.exons:
      if p < c + e.stop - e.start:
        # the position is within this exon
        return p - c + e.start
      else:
        # sorry mario, your position is in another exon
        c += e.stop - e.start
    assert(False)  # we should never get here

  def chromosomeCoordinateToExon(self, p):
    """ Take position P with 0-based chromosome-relative position and convert it
    to 0-based exon-relative position. If position does not exist in
    exon, return None.
    """
    if p is None: return None
    if self.chromosomeInterval.strand:
      def _stranded(v): return v
    else:
      def _stranded(v):
        return sum([(e.stop - e.start) for e in self.exons]) - 1 - v
    c = 0  # cumulative position through exon space
    e_start = self.exons[0].start
    for e in self.exons:
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

  def chromosomeCoordinateToMRna(self, p):
    """ Take position P with 0-based chromosome-relative position and convert it
    to 0-based mRNA-relative position. If position does not exist in mRNA,
    return None.
    """
    if p is None: return None
    if p < 0:
      return None
    if p >= self.chromosomeInterval.stop:
      return None
    q = self.chromosomeCoordinateToExon(p)
    if q is None:
      return None
    if q < 0:
      return None
    if q >= sum([(e.stop - e.start) for e in self.exons]):
      return None
    return self.exonCoordinateToMRna(q)

  def getMRna(self, sequence):
    """ Return the mRNA sequence for the transcript (based on the exons) using
    a SEQUENCE object as the source for dna sequence.
    The returned sequence is in the correct 5'-3' orientation (i.e. it has
    been reverse complemented if necessary).
    """
    assert(self.chromosomeInterval.chromosome == sequence.name)
    assert(self.chromosomeInterval.stop <= sequence.getLength())
    s = ''
    # chromosome    ttttTTTTTTTTTTTtttt  t: thin T: THICK
    # exon            eeeeee eeee eee
    # mrna              mmmm mmmm m
    for e in self.exons:
      if self.thickStart < e.start and e.stop < self.thickEnd:
        # squarely in the CDS
        s += sequence.sliceSequence(e.start, e.stop)
      elif (e.start <= self.thickStart and e.stop < self.thickEnd
            and self.thickStart < e.stop):
        # thickStart marks the start of the mRNA
        s += sequence.sliceSequence(self.thickStart, e.stop)
      elif e.start <= self.thickStart and self.thickEnd <= e.stop:
        # thickStart and thickEnd mark the whole mRNA
        s += sequence.sliceSequence(self.thickStart, self.thickEnd)
      elif (self.thickStart < e.start and self.thickEnd <= e.stop
            and e.start < self.thickEnd):
        # thickEnd marks the end of the mRNA
        s += sequence.sliceSequence(e.start, self.thickEnd)
    if not self.chromosomeInterval.strand:
      s = reverseComplement(s)
    return s

  def getMRnaCodon(i, sequence):
    """ Given a nonnegative number I, and a SEQUENCE object,
    return the three nucleotide codon from the mRNA.
    """
    return self.getMRna(sequence)[i * 3:i * 3 + 3]

  def bedString(self):
    """ Write a transcript object to the given file.
    """
    strandChar = '-'
    if self.chromosomeInterval.strand:
        strandChar = '+'
    return '\t'.join(
      [self.chromosomeInterval.chromosome,
       str(self.chromosomeInterval.start),
       str(self.chromosomeInterval.stop),
       self.name, str(self.score), strandChar,
       str(self.thickStart), str(self.thickEnd),
       self.itemRgb, str(len(self.exons)),
       ','.join([str(exon.stop - exon.start) for exon in self.exons]),
       ','.join([str(exon.start - self.chromosomeInterval.start)
                 for exon in self.exons])])

  def getIntrons(self):
    """Get a list of ChromosomeIntervals representing the introns for this
    transcript. The introns are in *+ strand of CHROMOSOME* ordering,
    not the order that they appear in the transcript!"""
    introns = []
    prevExon = None
    for exon in self.exons:
      if prevExon is not None:
        assert exon.start > prevExon.stop
        assert exon.strand == prevExon.strand
        intron = ChromosomeInterval(exon.chromosome,
                                    prevExon.stop,
                                    exon.start,
                                    exon.strand)
        introns.append(intron)
      prevExon = exon
    return introns

  def __cmp__(self, transcript):
    return cmp((self.chromosomeInterval, self.name),
               (transcript.chromosomeInterval, transcript.name))


def DirType(d):
  """ given a string path to a directory, D, verify it can be used.
  """
  d = os.path.abspath(d)
  if not os.path.exists(d):
    raise ArgumentTypeError('DirType:%s does not exist' % d)
  if not os.path.isdir(d):
    raise ArgumentTypeError('DirType:%s is not a directory' % d)
  if os.access(d, os.R_OK):
    return d
  else:
    raise ArgumentTypeError('DirType:%s is not a readable dir' % d)


def FileType(f):
  """ given a string path to a file, F, verify it can be used.
  """
  f = os.path.abspath(f)
  if not os.path.exists(f):
    raise ArgumentTypeError('FileType:%s does not exist' % f)
  if not os.path.isfile(f):
    raise ArgumentTypeError('FileType:%s is not a regular file' % f)
  if os.access(f, os.R_OK):
    return f
  else:
    raise ArgumentTypeError('FileType:%s is not a readable file' % f)


def boilerplateArguments(extraArguments=None, extraChecks=None):
  """ start from scratch and hand back a valid and initialized args object.
  """
  from argparse import ArgumentParser
  parser = ArgumentParser()
  initializeArguments(parser)
  if extraArguments is not None:
    extraArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  if extraChecks is not None:
    extraChecks(args, parser)
  return args


def initializeArguments(parser):
  """ given an argparse ArgumentParser object, add in the default arguments.
  """
  parser.add_argument('--refGenome', type=str)
  parser.add_argument('--genome', type=str)
  parser.add_argument('--geneCheckBed', type=FileType)
  parser.add_argument('--geneCheckBedDetails', type=FileType)
  parser.add_argument('--originalGeneCheckBed', type=FileType)
  parser.add_argument('--originalGeneCheckBedDetails', type=FileType)
  parser.add_argument('--alignment', type=FileType)
  parser.add_argument('--sequence', type=FileType)
  parser.add_argument('--refSequence', type=FileType)
  parser.add_argument('--chromSizes', type=FileType)
  parser.add_argument('--outDir', type=DirType)
  parser.add_argument('--mode', choices=['transmap', 'augustus'])


def checkArguments(args, parser):
  """ Make sure all of the args are properly set for the default arguments.
  """
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['refGenome', 'genome',
                 'geneCheckBed', 'geneCheckBedDetails',
                 'originalGeneCheckBed', 'originalGeneCheckBedDetails',
                 'alignment', 'sequence', 'refSequence', 'chromSizes',
                 'outDir'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)
  # record the issuing command
  with open(os.path.join(args.outDir, 'command.log'), 'w') as f:
    f.write('%s' % sys.argv[0])
    for name, value in pairs:
      f.write(' --%s %s' % (name, value))
    f.write('\n')


_nuc_pairs = [('a', 't'), ('g', 'c'), ('n', 'n')]
_complement = {}
for a, b in _nuc_pairs:
  _complement[a], _complement[b] = b, a
  _complement[a.upper()], _complement[b.upper()] = b.upper(), a.upper()
_complement['-'] = '-'
def complement(seq):
  """ given a sequence, return the complement.
  """
  seq = ''.join([_complement[s] for s in seq])
  return seq


def reverseComplement(seq):
  """ Given a sequence, return the reverse complement.
  """
  seq = seq[::-1]  # reverse
  seq = complement(seq)
  return seq


_codonToAminoAcid = {
  'ATG': 'Met',
  'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop', 'TAR': 'Stop', 'TRA': 'Stop',
  'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GCN': 'Ala',
  'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg','AGA': 'Arg',
  'AGG': 'Arg', 'CGN': 'Arg', 'MGR': 'Arg',
  'AAT': 'Asn', 'AAC': 'Asn', 'AAY': 'Asn',
  'GAT': 'Asp', 'GAC': 'Asp', 'GAY': 'Asp',
  'TGT': 'Cys', 'TGC': 'Cys', 'TGY': 'Cys',
  'CAA': 'Gln', 'CAG': 'Gln', 'CAR': 'Gln',
  'GAA': 'Glu', 'GAG': 'Glu', 'GAR': 'Glu',
  'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly', 'GGN': 'Gly',
  'CAT': 'His', 'CAC': 'His', 'CAY': 'His',
  'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATH': 'Ile',
  'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu',
  'CTG': 'Leu', 'YTR': 'Leu', 'CTN': 'Leu',
  'AAA': 'Lys', 'AAG': 'Lys', 'AAR': 'Lys',
  'TTT': 'Phe', 'TTC': 'Phe', 'TTY': 'Phe',
  'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CCN': 'Pro',
  'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser',
  'AGC': 'Ser', 'TCN': 'Ser', 'AGY': 'Ser',
  'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'ACN': 'Thr',
  'TGG': 'Trp',
  'TAT': 'Tyr', 'TAC': 'Tyr', 'TAY': 'Tyr',
  'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GTN': 'Val',
  }
def codonToAminoAcid(c):
  """ Given a codon C, return an amino acid or ??? if codon unrecognized.
  Codons could be unrecognized due to ambiguity IUPAC characters.
  """
  if c is None: return None
  c = c.upper()
  if c in _codonToAminoAcid:
    return _codonToAminoAcid[c]
  return '???'


def translateSequence(seq):
  """ Convert an entire DNA sequence to an amino acid sequence.
  """
  aa = ''
  for i in xrange(0, len(seq), 3):
    aa += codonToAminoAcid(seq[i:i+3])
  return aa


def readCodons(seq):
  """ Provide an iterator that reads through a sequence one codon at a time.
  """
  i = 0
  while i < len(seq):
    t = seq[i:i+3]
    i += 3
    yield t


def getSequences(infile, upper=False):
  """ Given a path to a fasta file, return a dictionary of Sequence objects
  keyed on the sequence name.
  """
  seqDict = {}
  seq = None
  with open(infile, 'r') as f:
    for seq in readSequence(f):
      seqDict[seq.name] = seq
      if upper:
        seqDict[seq.name].setUpper()
  return seqDict


def readSequence(infile):
  """ provide an iterator that reads through fasta files.
  """
  buff = None
  eat_buffer = True
  while True:
    if eat_buffer:
      # new record
      if buff is None:
        header = ''
        while not header.startswith('>'):
          header = infile.readline().strip()
          if header == '':
            return
      else:
        header = buff
      assert(header.startswith('>'))
      name = header.replace('>', '').strip().split(" ")[0]
      seq = ''
    line = infile.readline().strip()
    if line:
      if line.startswith('>'):
        # stop processing the record, store this line.
        buff = line
        eat_buffer = True
        yield Sequence(name, seq)
      else:
        eat_buffer = False
        seq += line
    else:
      # eof
      if buff is not None:
        buff = None
        yield Sequence(name, seq)
      else:
        if seq != '':
          yield Sequence(name, seq)
          name = ''
          seq = ''
        else:
          return


def getChromSizes(infile):
  """ read a chrom sizes file and return a dict keyed by names valued by ints.
  """
  chromDict = {}
  with open(infile, 'r') as f:
    for line in f:
      line = line.strip()
      if line == '':
        continue
      data = line.split()
      chromDict[data[0]] = int(data[1])
  return chromDict


def getAlignments(infile):
  """ read a PSL file and return a list of PslRow objects
  """
  psls = []
  with open(infile, 'r') as f:
    for psl in readPsls(f):
      psls.append(psl)
  return psls


def readPsls(infile):
  """ provide an iterator that reads through psl files.
  """
  while True:
    line = infile.readline().strip()
    if line == '':
      return
    yield PslRow(line)


def getTranscripts(bedFile, bedDetailsFile):
  """ Given a path to a standard BED file and a details BED, return a list of
  Transcript objects.
  """
  transcripts = []
  bf = open(bedFile, 'r')
  bdf = open(bedDetailsFile, 'r')
  for t in transcriptIterator(bf, bdf):
    transcripts.append(t)
  return transcripts


def transcriptListToDict(transcripts, noDuplicates=False):
  """ Given a list af Transcript objects, attempt to transfrom them into a dict
  of lists. key is transcript name, value is list of Transcript objects.
  If NODUPLICATES is true, then the value will be a single Transcript object.
  """
  result = {}
  for t in transcripts:
    if t.name not in result:
      result[t.name] = []
    else:
      if noDuplicates:
        raise RuntimeError('transcriptListToDict: Discovered a '
                           'duplicate transcript %s %s'
                           % (t.name, t.chromosomeInterval.chromosome))
    if noDuplicates:
      result[t.name] = t
    else:
      result[t.name].append(t)
  return result


def tokenizeBedStream(bedStream):
  """ Iterator through bed file, returning lines as list of tokens
  """
  for line in bedStream:
    if line != '':
      tokens = line.split()
      yield tokens


def normalizeAnnotation(transcriptAnnotation):
  """ Normalizes the transcript annotation labels.
  This is meant to munge the labels of transcript annotations according to
  hacky needs of the input annotation type labels.
  """
  if (len(transcriptAnnotation.labels) > 1 and
      (transcriptAnnotation.labels[0].endswith('Splice') and
       len(transcriptAnnotation.labels[1].split('.')) == 3)):
    # try to find lists like ['unknownUtrSplice', 'CC..AC']
    newLabels = [ '_'.join(transcriptAnnotation.labels[:2])]
    newLabels += transcriptAnnotation.labels[2:]
    transcriptAnnotation.labels = newLabels
  elif (len(transcriptAnnotation.labels) > 1 and
        transcriptAnnotation.labels[0] == 'orfStop' and
        transcriptAnnotation.labels[1] in ['TAA', 'TAG', 'TGA']):
    # try to find lists like ['orfStop', 'TAG']
    newLabels = [ '_'.join(transcriptAnnotation.labels[:2])]
    newLabels += transcriptAnnotation.labels[2:]
    transcriptAnnotation.labels = newLabels


def transcriptIterator(transcriptsBedStream, transcriptDetailsBedStream):
  """ Iterates over the transcripts detailed in the two streams, producing
  Transcript objects. Streams are any iterator that returns bedlines or empty
  strings.
  """
  transcriptsAnnotations = {}
  for tokens in tokenizeBedStream(transcriptDetailsBedStream):
    assert (len(tokens) == 4 or len(tokens) == 9)  # 9 if it has color data.
    tA = TranscriptAnnotation(
      ChromosomeInterval(tokens[0], tokens[1], tokens[2], None),
      tokens[3].split('/')[-1], tokens[3].split('/')[:-1])
    # normalizeAnnotation(tA)  # removed this to improve xml
    key = (tA.name, tA.chromosomeInterval.chromosome)
    if key not in transcriptsAnnotations:
      transcriptsAnnotations[key] = []
    transcriptsAnnotations[key].append(tA)

  for tokens in tokenizeBedStream(transcriptsBedStream):
    assert len(tokens) == 12
    # Transcript
    name = tokens[3]
    # Get the chromosome interval
    assert tokens[5] in ['+', '-']
    cI = ChromosomeInterval(tokens[0], tokens[1], tokens[2], tokens[5] == '+')
    # Get the exons
    def getExons(exonNumber, blockSizes, blockStarts):
      assert exonNumber == len(blockSizes)
      assert exonNumber == len(blockStarts)
      return [ChromosomeInterval(
          cI.chromosome, cI.start + int(blockStarts[i]),
          cI.start + int(blockStarts[i]) + int(blockSizes[i]), cI.strand)
              for i in range(exonNumber)]
    exons = getExons(int(tokens[9]),
                     tokens[10].split(','), tokens[11].split(','))
    # Get the name annotations
    annotations = []
    key = (name, cI.chromosome)
    if key in transcriptsAnnotations:
      annotations = transcriptsAnnotations[key]
    filteredAnnotations = []
    for tA in annotations:
      if cI.contains(tA.chromosomeInterval):
        tA.chromosomeInterval.strand = cI.strand
        filteredAnnotations.append(tA)
    yield Transcript(
      cI, name, exons, filteredAnnotations,
      int(tokens[4]), int(tokens[6]),
      int(tokens[7]), tokens[8])


def writeAllBeds(transcripts, args):
  """ Convenience function to take a list of TRANSCRIPTS and write out the
  standard and details beds to the expected location in args.outDir.
  """
  out_bed, out_bed_details = getBedOutFiles(args)
  writeTranscriptBedFile(transcripts, out_bed)
  writeDetailsBedFile(transcripts, out_bed_details)


def getBedOutFiles(args):
  """ Return the paths of the new bed and bedDetail files.
  """
  bed = os.path.join(args.outDir, 'out.bed')
  bedDetails = os.path.join(args.outDir, 'out_details.bed')
  return bed, bedDetails


def writeDetailsBedFile(transcripts, detailsBedFile):
  """ Writes out a details bed file for a set of transcripts - that is the set
  of annotations of the transcripts. The bed file must be in chromosome order.
  """
  annotations = list(
    set(reduce(lambda x, y : x + y,
               [transcript.annotations for transcript in transcripts])))
  annotations.sort()
  annotationsFileHandle = open(detailsBedFile, 'w')
  for annotation in annotations:
    annotationsFileHandle.write(annotation.bedString() + '\n')
  annotationsFileHandle.close()


def writeTranscriptBedFile(transcripts, bedFile):
  """ Writes out an bed file for a set of transcripts.
  """
  transcripts = transcripts[:]  # no side effects
  transcripts.sort()
  bedFileHandle = open(bedFile, 'w')
  for transcript in transcripts:
    bedFileHandle.write(transcript.bedString() + '\n')
  bedFileHandle.close()


def removeAlignmentNumber(s):
  """ If the name of the transcript ends with -d as in
  ENSMUST00000169901.2-1, return ENSMUST00000169901.2
  """
  s = s[:]
  i = s.find('-')
  if i == -1:
    return s
  else:
    return s[0:i]
