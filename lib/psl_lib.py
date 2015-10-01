"""
Convenience library for working with psl alignment files.

Original Author: Dent Earl
Modified by Ian Fiddes
"""

from collections import Counter
import re

__author__ = "Ian Fiddes"


class PslRow(object):
    """ Represents a single row in a PSL file.
    http://genome.ucsc.edu/FAQ/FAQformat.html#format2
    """
    __slots__ = ('matches', 'mismatches', 'repmatches', 'n_count', 'q_num_insert', 'q_base_insert', 't_num_insert',
                 't_base_insert', 'strand', 'q_name', 'q_size', 'q_start', 'q_end', 't_name', 't_size', 't_start',
                 't_end', 'block_count', 'block_sizes', 'q_starts', 't_starts')

    def __init__(self, line):
        data = line.split()
        assert(len(data) == 21)
        self.matches = int(data[0])
        self.mismatches = int(data[1])
        self.repmatches = int(data[2])
        self.n_count = int(data[3])
        self.q_num_insert = int(data[4])
        self.q_base_insert = int(data[5])
        self.t_num_insert = int(data[6])
        self.t_base_insert = int(data[7])
        self.strand = data[8]
        self.q_name = data[9]
        self.q_size = int(data[10])
        self.q_start = int(data[11])
        self.q_end = int(data[12])
        self.t_name = data[13]
        self.t_size = int(data[14])
        self.t_start = int(data[15])
        self.t_end = int(data[16])
        self.block_count = int(data[17])
        # lists of ints
        self.block_sizes = [int(x) for x in data[18].split(',') if x]
        self.q_starts = [int(x) for x in data[19].split(',') if x]
        self.t_starts = [int(x) for x in data[20].split(',') if x]

    def hash_key(self):
        """ return a string to use as dict key.
        """
        return '%s_%s_%d_%d' % (self.q_name, self.t_name, self.t_start, self.t_end)

    def target_coordinate_to_query(self, p):
        """ Take position P in target coordinates (positive) and convert it
        to query coordinates (positive).
        """
        if p < self.t_start:
            return None
        if p >= self.t_end:
            return None
        if self.strand not in ['+', '-']:
            raise RuntimeError('Unanticipated strand: %s' % self.strand)
        for i, t in enumerate(self.t_starts):
            if p < t:
                continue
            if p >= t + self.block_sizes[i]:
                continue
            # p must be in block
            offset = p - t
            if self.strand == '+':
                return self.q_starts[i] + offset
            else:
                return self.q_size - (self.q_starts[i] + offset) - 1
        return None

    def query_coordinate_to_target(self, p):
        """ Take position P in query coordinates (positive) and convert it
        to target coordinates (positive).
        """
        if p < self.q_start:
            return None
        if p >= self.q_end:
            return None
        if self.strand not in ['+', '-']:
            raise RuntimeError('Unanticipated strand: %s' % self.strand)
        # this is the easier one to write
        if self.strand == '-':
            p = self.q_size - p - 1
        for i, q in enumerate(self.q_starts):
            if p < q:
                continue
            if p >= q + self.block_sizes[i]:
                continue
            # p must be in block
            offset = p - q
            return self.t_starts[i] + offset
        return None

    def psl_string(self):
        """ return SELF as a psl formatted line.
        """
        return "\t".join([self.matches, self.mismatches, self.repmatches, self.n_count, self.q_num_insert,
                          self.q_base_insert, self.t_num_insert, self.t_base_insert, self.strand, self.q_name,
                          self.q_start, self.q_end, self.t_name, self.t_size, self.t_start, self.t_end,
                          self.block_count, ','.join([str(b) for b in self.block_sizes]),
                          ','.join([str(b) for b in self.q_starts]),
                          ','.join([str(b) for b in self.t_starts])])


def read_psl(infile, uniqify=False):
    """ read a PSL file and return a list of PslRow objects
    """
    psls = []
    with open(infile, 'r') as f:
        for psl in psl_iterator(f, uniqify):
            psls.append(psl)
    return psls


def psl_iterator(infile, uniqify=False):
    """ Iterator to loop over psls returning PslRow objects.
    If uniqify is set, will add a number to each name starting with -1"""
    names = Counter()
    while True:
        line = infile.readline().strip()
        if line == '':
            return
        r = PslRow(line)
        if uniqify is False:
            yield r
        else:
            names[remove_alignment_number(r.q_name)] += 1
            yield uniqify_psl_row(r, names[r.q_name])


def get_psl_dict(alignments):
    """
    turns an alignment list from readPsl to a dict keyed on alignmentID.
    """
    alignments_dict = {}
    for a in alignments:
        if a.qName in alignments_dict:
                raise RuntimeError("getPslDict found duplicate transcript {} when noDuplicates was set".format(a.qName))
        else:
            alignments_dict[a.qName] = a
    return alignments_dict


def remove_alignment_number(s, aln_re=re.compile("-[0-9]+$")):
    """ If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    """
    return aln_re.split(s)[0]


# aug_re = re.compile("^((aug-[0-9]+)|(aug))-")
# new regular expression to match new naming scheme for 1509
def remove_augustus_alignment_number(s, aug_re=re.compile("^((augI[0-9]+-[0-9]+)|(augI[0-9]+))-")):
    """
    removes the alignment numbers prepended by augustus
    """
    return aug_re.split(s)[-1]


def uniqify_psl_row(row, val):
    """ Uniqifies the name of <row> by adding -<val> to it
    """
    row.qName = "-".join([row.qName, str(val)])
    return row