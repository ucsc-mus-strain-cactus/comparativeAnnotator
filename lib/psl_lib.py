"""
Convenience library for working with psl alignment files.

Original Author: Dent Earl
Modified by Ian Fiddes
"""


from collections import defaultdict, Counter

class PslRow(object):
    """ Represents a single row in a PSL file.
    http://genome.ucsc.edu/FAQ/FAQformat.html#format2
    """
    __slots__ = ('matches', 'misMatches', 'repMatches', 'nCount',
                     'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert',
                     'strand', 'qName', 'qSize', 'qStart', 'qEnd',
                     'tName', 'tSize', 'tStart', 'tEnd', 'blockCount',
                     'blockSizes', 'qStarts', 'tStarts') # conserve memory
    
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


def readPsl(infile, uniqify=False):
    """ read a PSL file and return a list of PslRow objects
    """
    psls = []
    with open(infile, 'r') as f:
        for psl in pslIterator(f, uniqify):
            psls.append(psl)
    return psls


def pslIterator(infile, uniqify=False):
    """ Iterator to loop over psls returning PslRow objects.
    If uniqify is set, will add a number to each name starting with -0"""
    names = Counter()
    while True:
        line = infile.readline().strip()
        if line == '':
            return
        r = PslRow(line)
        if uniqify is False:
            yield r
        else:
            yield uniqifyPslRow(r, names[r.qName])
            names[removeAlignmentNumber(r.qName)] += 1


def getPslDict(alignments, noDuplicates=False):
    """
    turns an alignment list from readPsl to a dict keyed on alignmentID.
    """
    alignments_dict = {}
    for a in alignments:
        if a.qName not in alignments_dict:
            alignments_dict[a.qName] = []
        else:
            if noDuplicates:
                raise RuntimeError("getPslDict found duplicate transcript {}"
                    "when noDuplicates was set".format(a.qName))
        if noDuplicates:
            alignments_dict[a.qName] = a
        else:
            alignments_dict[a.qName].append(a)
    return alignments_dict


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


def uniqifyPslRow(row, val):
    """ Uniqifies the name of <row> by adding -<val> to it
    """
    row.qName = "-".join([row.qName, str(val)])
    return row