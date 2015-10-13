"""
Convenience library for sequence information, including BED/genePred files and fasta files
Needs the python library pyfaidx installed.

Original Author: Dent Earl
Modified by Ian Fiddes
"""

import string
import copy
import math
from itertools import izip
from lib.general_lib import tokenize_stream
from pyfaidx import Fasta

__author__ = "Ian Fiddes"


class Transcript(object):
    """
    Represent a transcript record from a bed file. Stores the fields from the BED file
    and then uses them to create the following class members:
    chromosomeInterval: a ChromosomeInterval object representing the entire transcript
        in chromosome coordinates.
    exonIntervals: a list of ChromosomeInterval objects representing each exon in
        chromosome coordinates.
    intronIntervals: a list of ChromosomeInterval objects representing each intron
        in chromosome coordinates.
    exons: a list of Exon objects representing this transcript. These objects store mappings
        between chromosome, transcript and CDS coordinate space. Transcript and CDS coordinates
        are always transcript relative (5'->3').

    To be more efficient, the cds and mRNA slots are saved for if those sequences are ever retrieved.
    Then they will be stored so we don't slice the same thing over and over.
    """

    __slots__ = ('name', 'strand', 'score', 'thick_start', 'rgb', 'thick_stop', 'start', 'stop', 'intron_intervals',
                 'exon_intervals', 'exons', 'cds', 'mrna', 'block_sizes', 'block_starts', 'block_count', 'chromosome',
                 'cds_size', 'transcript_size')

    def __init__(self, bed_tokens):
        self.chromosome = bed_tokens[0]
        self.start = int(bed_tokens[1])
        self.stop = int(bed_tokens[2])
        self.name = bed_tokens[3]
        self.score = int(bed_tokens[4])
        self.strand = convert_strand(bed_tokens[5])
        self.thick_start = int(bed_tokens[6])
        self.thick_stop = int(bed_tokens[7])
        self.rgb = bed_tokens[8]
        self.block_count = bed_tokens[9]
        self.block_sizes = bed_tokens[10]
        self.block_starts = bed_tokens[11]
        # build chromosome intervals for exons and introns
        self.exon_intervals = self._get_exon_intervals(bed_tokens)
        self.intron_intervals = self._get_intron_intervals()
        # build Exons mapping transcript space coordinates to chromosome
        self.exons = self._get_exons(bed_tokens)
        # calculate sizes
        self._get_cds_size()
        self._get_size()

    def __len__(self):
        return self.transcript_size

    def __cmp__(self, transcript):
        return cmp((self, self.name), (transcript, transcript.name))

    def get_bed(self, rgb=None, name=None, start_offset=None, stop_offset=None):
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
            return [self.chromosome, self.start, self.stop, name, self.score, convert_strand(self.strand),
                    self.thick_start, self.thick_stop, rgb, self.block_count, self.block_sizes, self.block_starts]
        elif start_offset == stop_offset:
            assert self.chromosome_coordinate_to_transcript(start_offset) is not None   # no intron records
            return [self.chromosome, start_offset, stop_offset, name, self.score, convert_strand(self.strand),
                    start_offset, stop_offset, rgb, 1, 0, 0]

        def _move_start(exon_intervals, block_count, block_starts, block_sizes, start, start_offset):
            to_remove = len([x for x in exon_intervals if x.start <= start_offset and x.stop <= start_offset])
            assert to_remove < len(exon_intervals)
            if to_remove > 0:
                block_count -= to_remove
                block_sizes = block_sizes[to_remove:]
                start += block_starts[to_remove]
                new_block_starts = [0]
                for i in xrange(to_remove, len(block_starts) - 1):
                    new_block_starts.append(block_starts[i + 1] - block_starts[i] + new_block_starts[-1])
                block_starts = new_block_starts
            if start_offset > start:
                block_sizes[0] += start - start_offset
                block_starts[1:] = [x + start - start_offset for x in block_starts[1:]]
                start = start_offset
            return start, block_count, block_starts, block_sizes

        def _move_stop(exon_intervals, block_count, block_starts, block_sizes, stop, start, stop_offset):
            to_remove = len([x for x in exon_intervals if x.stop >= stop_offset and x.start >= stop_offset])
            assert to_remove < len(exon_intervals)
            if to_remove > 0:
                block_count -= to_remove
                block_sizes = block_sizes[:-to_remove]
                block_starts = block_starts[:-to_remove]
                assert len(block_sizes) == len(block_starts)
                if len(block_sizes) == 0:
                    block_sizes = block_starts = [0]
                    block_count = 1
                stop = start + block_sizes[-1] + block_starts[-1]
            if start + block_starts[-1] < stop_offset < stop:
                block_sizes[-1] = stop_offset - start - block_starts[-1]
                stop = stop_offset
            return stop, block_count, block_starts, block_sizes

        block_count = int(self.block_count)
        block_starts = map(int, self.block_starts.split(","))
        block_sizes = map(int, self.block_sizes.split(","))
        start = self.start
        stop = self.stop
        thick_start = self.thick_start
        thick_stop = self.thick_stop

        if start_offset is not None and start_offset > start:
            start, block_count, block_starts, block_sizes = _move_start(self.exon_intervals, block_count, block_starts,
                                                                    block_sizes, start, start_offset)
        if stop_offset is not None and stop_offset < stop:
            stop, block_count, block_starts, block_sizes = _move_stop(self.exon_intervals, block_count, block_starts,
                                                                  block_sizes, stop, start, stop_offset)
        if start > thick_start:
            thick_start = start
        if stop < thick_stop:
            thick_stop = stop
        if (start > thick_stop and stop > thick_stop) or (start < thick_start and stop < thick_start):
            thick_start = 0
            thick_stop = 0
        block_starts = ",".join(map(str, block_starts))
        block_sizes = ",".join(map(str, block_sizes))
        return [self.chromosome, start, stop, name, self.score, convert_strand(self.strand), thick_start, thick_stop, rgb,
                block_count, block_sizes, block_starts]

    def _get_exon_intervals(self, bed_tokens):
        """
        Gets a list of exon intervals in chromosome coordinate space.
        These exons are on (+) strand ordering regardless of transcript strand.
        This means (-) strand genes will be represented backwards
        """
        exons = []
        start, stop = int(bed_tokens[1]), int(bed_tokens[2])
        chrom, strand = bed_tokens[0], convert_strand(bed_tokens[5])

        block_sizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        block_starts = [int(x) for x in bed_tokens[11].split(",") if x != ""]

        for block_size, block_start in izip(block_sizes, block_starts):
            exons.append(ChromosomeInterval(chrom, start + block_start, start + block_start + block_size, strand))
        return exons

    def _get_intron_intervals(self):
        """
        Get a list of ChromosomeIntervals representing the introns for this
        transcript. The introns are in *+ strand of CHROMOSOME* ordering,
        not the order that they appear in the transcript!
        """
        introns = []
        prev_exon = None
        for exon in self.exon_intervals:
            if prev_exon is not None:
                assert exon.strand == prev_exon.strand
                assert exon.start >= prev_exon.stop
                intron = ChromosomeInterval(exon.chromosome, prev_exon.stop, exon.start, exon.strand)
                introns.append(intron)
            prev_exon = exon
        return introns

    def _get_exons(self, bed_tokens):
        """
        Get a list of Exons representing the exons in transcript coordinate
        space. This is in transcript order. See the Exon class for more.
        """
        exons = []
        chrom_start, chrom_stop = int(bed_tokens[1]), int(bed_tokens[2])
        thick_start, thick_stop = int(bed_tokens[6]), int(bed_tokens[7])
        if thick_start == thick_stop:
            thick_start = thick_stop = 0
        chrom, strand = bed_tokens[0], convert_strand(bed_tokens[5])

        block_count = int(bed_tokens[9])
        block_sizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        block_starts = [int(x) for x in bed_tokens[11].split(",") if x != ""]

        ##################################################################
        # HERE BE DRAGONS
        # this is seriously ugly code to maintain proper mapping
        # between coordinate spaces. See the unit tests.
        ##################################################################
        if strand is False:
            block_sizes = reversed(block_sizes)
            block_starts = reversed(block_starts)

        t_pos, cds_pos = 0, None
        for block_size, block_start in izip(block_sizes, block_starts):
            # calculate transcript relative coordinates
            this_start = t_pos
            this_stop = t_pos + block_size
            # calculate chromosome relative coordinates
            this_chrom_start = chrom_start + block_start
            this_chrom_stop = chrom_start + block_start + block_size
            # calculate transcript-relative CDS positions
            # cds_pos is pos of first coding base in CDS coordinates
            this_cds_start, this_cds_stop, this_cds_pos = None, None, None
            if strand is True:
                # special case - single exon
                if block_count == 1:
                    this_cds_pos = 0
                    this_cds_start = thick_start - this_chrom_start
                    this_cds_stop = thick_stop - this_chrom_start
                # special case - entirely non-coding
                elif thick_start == thick_stop == 0:
                    this_cds_start, this_cds_stop, this_cds_pos = None, None, None
                # special case - CDS starts and stops on the same exon
                elif (this_chrom_start <= thick_start < this_chrom_stop and this_chrom_start < thick_stop <=
                        this_chrom_stop):
                    this_cds_pos = 0
                    cds_pos = this_chrom_stop - thick_start
                    this_cds_start = this_start + thick_start - this_chrom_start
                    this_cds_stop = this_stop + thick_stop - this_chrom_stop
                # is this the start codon containing exon?
                elif this_chrom_start <= thick_start < this_chrom_stop:
                    cds_pos = this_chrom_stop - thick_start
                    this_cds_pos = 0
                    this_cds_start = this_start + thick_start - this_chrom_start
                # is this the stop codon containing exon?
                elif this_chrom_start < thick_stop <= this_chrom_stop:
                    this_cds_pos = cds_pos
                    cds_pos += thick_stop - this_chrom_start
                    this_cds_stop = this_stop + thick_stop - this_chrom_stop
                # is this exon all coding?
                elif (this_cds_stop is None and this_cds_start is None and thick_stop >=
                      this_chrom_stop and thick_start < this_chrom_start):
                    this_cds_pos = cds_pos
                    cds_pos += block_size
            else:
                # special case - single exon
                if block_count == 1:
                    this_cds_pos = 0
                    this_cds_start = this_chrom_stop - thick_stop
                    this_cds_stop = thick_stop - this_chrom_start + this_cds_start
                # special case - entirely non-coding
                elif thick_start == thick_stop == 0:
                    this_cds_start, this_cds_stop, this_cds_pos = None, None, None
                # special case - start and stop codons are on the same exon
                elif (this_chrom_start < thick_stop <= this_chrom_stop and this_chrom_start <= thick_start <
                      this_chrom_stop):
                    cds_pos = thick_stop - this_chrom_start
                    this_cds_pos = 0
                    this_cds_start = this_start + this_chrom_stop - thick_stop
                    this_cds_stop = this_start + this_chrom_stop - thick_start
                # is this the start codon containing exon?
                elif this_chrom_start < thick_stop <= this_chrom_stop:
                    cds_pos = thick_stop - this_chrom_start
                    this_cds_pos = 0
                    this_cds_start = this_start + this_chrom_stop - thick_stop
                # is this the stop codon containing exon?
                elif this_chrom_start <= thick_start < this_chrom_stop:
                    this_cds_pos = cds_pos
                    this_cds_stop = this_start + this_chrom_stop - thick_start
                # is this exon all coding?
                elif (this_cds_stop is None and this_cds_start is None and thick_stop >=
                      this_chrom_stop and thick_start < this_chrom_start):
                    this_cds_pos = cds_pos
                    cds_pos += block_size
            exons.append(Exon(this_start, this_stop, strand, this_chrom_start, this_chrom_stop, this_cds_start,
                              this_cds_stop, this_cds_pos))
            t_pos += block_size
        return exons

    def _get_size(self):
        self.transcript_size = sum(x.stop - x.start for x in self.exon_intervals)

    def _get_cds_size(self):
        l = 0
        for e in self.exon_intervals:
            if self.thick_start < e.start and e.stop < self.thick_stop:
                # squarely in the CDS
                l += e.stop - e.start
            elif e.start <= self.thick_start < e.stop < self.thick_stop:
                # thickStart marks the start of the CDS
                l += e.stop - self.thick_start
            elif e.start <= self.thick_start and self.thick_stop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                l += self.thick_stop - self.thick_start
            elif self.thick_start < e.start < self.thick_stop <= e.stop:
                # thickStop marks the end of the CDS
                l += self.thick_stop - e.start
        self.cds_size = l

    def get_mrna(self, seq_dict):
        """
        Returns the mRNA sequence for this transcript based on a Fasta object.
        and the start/end positions and the exons. Sequence returned in
        5'-3' transcript orientation.
        """
        if hasattr(self, "mrna"):
            return self.mrna
        sequence = seq_dict[self.chromosome]
        assert self.stop <= len(sequence)
        s = []
        for e in self.exon_intervals:
            s.append(sequence[e.start:e.stop])
        if self.strand is True:
            mrna = "".join(s)
        else:
            mrna = reverse_complement("".join(s))
        self.mrna = mrna.upper()
        return mrna

    def get_sequence(self, seq_dict):
        """
        Returns the entire chromosome sequence for this transcript, (+) strand orientation.
        """
        sequence = seq_dict[self.chromosome]
        return sequence[self.start:self.stop].upper()

    def get_cds(self, seq_dict):
        """
        Return the CDS sequence (as a string) for the transcript
        (based on the exons) using a sequenceDict as the sequence source.
        The returned sequence is in the correct 5'-3' orientation (i.e. it has
        been reverse complemented if necessary).
        """
        if hasattr(self, "cds"):
            return self.cds
        sequence = seq_dict[self.chromosome]
        assert self.stop <= len(sequence)
        # make sure this isn't a non-coding gene
        if self.thick_start == self.thick_stop == 0:
            return ""
        s = []
        for e in self.exon_intervals:
            if self.thick_start < e.start and e.stop < self.thick_stop:
                # squarely in the CDS
                s.append(sequence[e.start:e.stop])
            elif e.start <= self.thick_start < e.stop < self.thick_stop:
                # thickStart marks the start of the CDS
                s.append(sequence[self.thick_start:e.stop])
            elif e.start <= self.thick_start and self.thick_stop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                s.append(sequence[self.thick_start: self.thick_stop])
            elif self.thick_start < e.start < self.thick_stop <= e.stop:
                # thickStop marks the end of the CDS
                s.append(sequence[e.start:self.thick_stop])
        if not self.strand:
            cds = reverse_complement("".join(s)).upper()
        else:
            cds = "".join(s)
        self.cds = cds.upper()
        return cds

    def get_transcript_coordinate_cds_start(self):
        """
        Returns the transcript-relative position of the CDS start
        """
        return self.chromosome_coordinate_to_transcript(self.thick_start)

    def get_transcript_coordinate_cds_stop(self):
        """
        Returns the transcript-relative position of the CDS stop
        """
        return self.chromosome_coordinate_to_transcript(self.thick_stop)

    def get_chromosome_coordinate_cds_start(self):
        """
        Returns the chromosome-relative position of the CDS start.
        This is not thickStart if a negative strand gene.
        Therefore, no guarantee it is smaller than the stop.
        """
        if self.strand is True:
            return self.thick_start
        else:
            return self.thick_stop - 1

    def get_chromosome_coordinate_cds_stop(self):
        """
        Returns the chromosome-relative position of the CDS stop.
        This is not thickStart if a negative strand gene.
        Therefore, no guarantee it is larger than the stop.
        """
        if self.strand is True:
            return self.thick_stop
        else:
            return self.thick_start - 1

    def get_protein_sequence(self, seq_dict):
        """
        Returns the translated protein sequence for this transcript in single
        character space.
        """
        cds = self.get_cds(seq_dict)
        if len(cds) < 3:
            return ""
        return translate_sequence(self.get_cds(seq_dict))

    def intron_sequence_iterator(self, seq_dict):
        """
        Iterates over intron sequences in transcript order and strand
        """
        if self.strand is True:
            for intron in self.intron_intervals:
                yield intron.get_sequence(seq_dict)
        else:
            for intron in reversed(self.intron_intervals):
                yield intron.get_sequence(seq_dict)

    def get_intron_sequences(self, seq_dict):
        """
        Wrapper for intronSequenceIterator that returns a list of sequences.
        """
        return list(self.intron_sequence_iterator(seq_dict))

    def transcript_coordinate_to_cds(self, p):
        """
        Takes a transcript-relative position and converts it to CDS coordinates.
        Will return None if this transcript coordinate is non-coding.
        Transcript/CDS coordinates are 0-based half open on 5'->3' transcript orientation.
        """
        for exon in self.exons:
            t = exon.transcript_pos_to_cds_pos(p)
            if t is not None:
                return t
        return None

    def transcript_coordinate_to_chromosome(self, p):
        """
        Takes a mRNA-relative position and converts it to chromosome position.
        Take a look at the docstring in the Exon class method chromPosToTranscriptPos
        for details on how this works.
        """
        for exon in self.exons:
            t = exon.transcript_pos_to_chrom_pos(p)
            if t is not None:
                return t
        return None

    def chromosome_coordinate_to_transcript(self, p):
        """
        Takes a chromosome-relative position and converts it to transcript
        coordinates. Transcript coordinates are 0-based half open on
        5'->3' transcript orientation.
        """
        for exon in self.exons:
            t = exon.chrom_pos_to_transcript_pos(p)
            if t is not None:
                return t
        return None

    def chromosome_coordinate_to_cds(self, p):
        """
        Takes a chromosome-relative position and converts it to CDS coordinates.
        Will return None if this chromosome coordinate is not in the CDS.
        """
        for exon in self.exons:
            t = exon.chrom_pos_to_cds_pos(p)
            if t is not None:
                return t
        return None

    def cds_coordinate_to_transcript(self, p):
        """
        Takes a CDS-relative position and converts it to Transcript coordinates.
        """
        for exon in self.exons:
            t = exon.cds_pos_to_transcript_pos(p)
            if t is not None:
                return t
        return None

    def cds_coordinate_to_chromosome(self, p):
        """
        Takes a CDS-relative position and converts it to Chromosome coordinates.
        """
        for exon in self.exons:
            t = exon.cds_pos_to_chrom_pos(p)
            if t is not None:
                return t
        return None

    def cds_coordinate_to_amino_acid(self, p, seq_dict):
        """
        Takes a CDS-relative position and a Fasta object that contains this
        transcript and returns the amino acid at that CDS position.
        Returns None if this is invalid.
        """
        cds = self.get_cds(seq_dict)
        if p >= len(cds) or p < 0:
            return None
        # we add 0.1 to the math.ceiling to make multiples of 3 work
        start, stop = int(math.floor(p / 3.0) * 3),  int(math.ceil((p + 0.1) / 3.0) * 3)
        if stop - start != 3:
            return None
        codon = cds[start: stop]
        return codon_to_amino_acid(codon)

    def transcript_coordinate_to_amino_acid(self, p, seq_dict):
        """
        Takes a transcript coordinate position and a Fasta object that contains
        this transcript and returns the amino acid at that transcript position.
        If this position is not inside the CDS returns None.
        """
        cds_pos = self.transcript_coordinate_to_cds(p)
        if cds_pos is None:
            return None
        return self.cds_coordinate_to_amino_acid(cds_pos, seq_dict)

    def chromosome_coordinate_to_amino_acid(self, p, seq_dict):
        """
        Takes a chromosome coordinate and a Fasta object that contains this
        transcript and returns the amino acid at that chromosome position.
        Returns None if this position is not inside the CDS.
        """
        cds_pos = self.chromosome_coordinate_to_cds(p)
        if cds_pos is None:
            return None
        return self.cds_coordinate_to_amino_acid(cds_pos, seq_dict)


class GenePredTranscript(Transcript):
    """
    Represent a transcript record from a genePred file. Stores the fields from the file
    and then uses them to create the following class members:
    chromosomeInterval: a ChromosomeInterval object representing the entire transcript
        in chromosome coordinates.
    exon_intervals: a list of ChromosomeInterval objects representing each exon in
        chromosome coordinates.
    intron_intervals: a list of ChromosomeInterval objects representing each intron
        in chromosome coordinates.
    exons: a list of Exon objects representing this transcript. These objects store mappings
        between chromosome, transcript and CDS coordinate space. Transcript and CDS coordinates
        are always transcript relative (5'->3').

    To be more efficient, the cds and mRNA slots are saved for if those sequences are ever retrieved.
    Then they will be stored so we don't slice the same thing over and over.
    """
    # adding slots for new fields
    __slots__ = ('cds_start_stat', 'cds_end_stat', 'exon_frames')

    def __init__(self, gene_pred_tokens):
        # Text genePred fields
        self.name = gene_pred_tokens[0]
        self.chromosome = gene_pred_tokens[1]
        self.strand = convert_strand(gene_pred_tokens[2])
        # Integer genePred fields
        self.score = 0  # no score in genePred files
        self.thick_start = int(gene_pred_tokens[5])
        self.thick_stop = int(gene_pred_tokens[6])
        self.start = int(gene_pred_tokens[3])
        self.stop = int(gene_pred_tokens[4])
        self.rgb = "128,0,0"  # no RGB in genePred files
        # genePred specific fields
        self.cds_start_stat = gene_pred_tokens[12]
        self.cds_end_stat = gene_pred_tokens[13]
        self.exon_frames = [int(x) for x in gene_pred_tokens[14].split(",") if x != ""]
        # convert genePred format coordinates to BED-like coordinates to make intervals
        self.block_count = gene_pred_tokens[7]
        block_starts = [int(x) for x in gene_pred_tokens[8].split(",") if x != ""]
        block_ends = [int(x) for x in gene_pred_tokens[9].split(",") if x != ""]
        self.block_sizes = ",".join(map(str, [e - s for e, s in izip(block_ends, block_starts)]))
        self.block_starts = ",".join(map(str, [x - self.start for x in block_starts]))
        bed_tokens = [gene_pred_tokens[1], self.start, self.stop, self.name, self.score, gene_pred_tokens[2],
                      self.thick_start, self.thick_stop, self.rgb, self.block_count,
                      self.block_sizes, self.block_starts]
        # build chromosome intervals for exons and introns
        self.exon_intervals = self._get_exon_intervals(bed_tokens)
        self.intron_intervals = self._get_intron_intervals()
        # build Exons mapping transcript space coordinates to chromosome
        self.exons = self._get_exons(bed_tokens)
        # calculate sizes
        self._get_cds_size()
        self._get_size()

    def get_protein_sequence(self, seq_dict):
        """
        Returns the translated protein sequence for this transcript in single
        character space. Overrides this function in the Transcript class to make use of frame information.
        """
        offset = find_offset(self.exon_frames, self.strand)
        cds = self.get_cds(seq_dict)
        if len(cds) < 3:
            return ""
        return translate_sequence(cds[offset:])


class Exon(object):
    """
    An Exon object stores information about one exon in both
    transcript coordinates (5'->3') and chromosome coordinates (+) strand.

    Transcript coordinates and chromosome coordinates are 0-based half open.

    Has methods to convert between the two coordinates, taking strand
    into account.

    If this exon contains the start or stop codon, contains those positions
    in transcript coordinates. If this exon is coding, contains the
    CDS-coordinate position of the first base.
    """
    __slots__ = ('start', 'stop', 'strand', 'chrom_start', 'chrom_stop', 'cds_start', 'cds_stop', 'cds_pos')

    def __init__(self, start, stop, strand, chrom_start, chrom_stop, cds_start, cds_stop, cds_pos):
        assert chrom_stop - chrom_start == stop - start
        self.strand = strand
        # start, stop are transcript-relative coordinates
        self.start = start
        self.stop = stop
        self.strand = strand
        self.chrom_start = chrom_start
        self.chrom_stop = chrom_stop
        # cdsStart/cdsStop are transcript-coordinate
        # None if not stop/start codon in this exon
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        # cdsPos is cds coordinate
        self.cds_pos = cds_pos

    def __len__(self):
        return self.stop - self.start

    def contains_chrom_pos(self, p):
        """does this exon contain a given chromosome position?"""
        if p is None:
            return None
        elif self.chrom_start <= p < self.chrom_stop:
            return True
        return False

    def contains_transcript_pos(self, p):
        """does this exon contain a given transcript position?"""
        if p is None:
            return None
        elif self.start <= p < self.stop:
            return True
        return False

    def contains_cds_pos(self, p):
        """does this exon contain a given CDS position?"""
        if p is None:
            return None
        # see cdsPosToTranscriptPos - it will return None on invalid CDS positions
        elif self.cds_pos_to_transcript_pos(p) is not None:
            return True
        return False

    def contains_cds(self):
        """does this exon contain CDS?"""
        if self.cds_start == self.cds_stop == self.cds_pos is None:
            return False
        return True

    def chrom_pos_to_transcript_pos(self, p):
        """
        Given a chromosome position, returns the transcript position
        if it is within this exon. Otherwise, returns None
        Chromosome position is always on (+) strand
        """
        if p is None:
            return None
        elif p < self.chrom_start or p >= self.chrom_stop:
            return None
        elif self.strand is True:
            return self.start + p - self.chrom_start
        else:
            return self.start + self.chrom_stop - 1 - p

    def chrom_pos_to_cds_pos(self, p):
        """
        Given a chromosome position, returns the CDS position if the given
        chromosome position is in fact a CDS position on this exon.
        """
        t_pos = self.chrom_pos_to_transcript_pos(p)
        return self.transcript_pos_to_cds_pos(t_pos)

    def transcript_pos_to_cds_pos(self, p):
        """
        Given a transcript position, report cds-relative position.
        Returns None if this transcript position is not coding.
        """
        if p is None:
            return None
        elif p < self.start or p >= self.stop:
            return None
        # is this a coding exon?
        elif self.contains_cds() is False:
            return None
        # special case of single exon gene
        if self.cds_start is not None and self.cds_stop is not None:
            if p < self.cds_start or p >= self.cds_stop:
                return None
            return p - self.cds_start
        # exon contains start codon
        elif self.cds_start is not None:
            # make sure we are within CDS
            if self.cds_start > p:
                return None
            return p - self.cds_start
        # exon contains stop codon
        elif self.cds_stop is not None:
            # make sure we are within CDS
            if self.cds_stop <= p:
                return None
            return self.cds_pos + p - self.start
        # exon must be entirely coding
        else:
            return self.cds_pos + p - self.start

    def transcript_pos_to_chrom_pos(self, p):
        """
        Given a transcript position, returns the chromosome position
        if it is within this exon otherwise return None
        """
        if p is None:
            return None
        elif p < self.start or p >= self.stop:
            return None
        elif self.strand is True:
            return p + self.chrom_start - self.start
        else:
            return self.chrom_stop + self.start - 1 - p

    def cds_pos_to_chrom_pos(self, p):
        """
        Given a cds position, returns the chromosome position
        if the cds position is on this exon
        """
        t_pos = self.cds_pos_to_transcript_pos(p)
        return self.transcript_pos_to_chrom_pos(t_pos)

    def cds_pos_to_transcript_pos(self, p):
        """
        Given a CDS position, returns the transcript position if it exists.
        Otherwise returns None.
        """
        if p is None:
            return None
        # not a coding exon
        elif self.contains_cds() is False:
            return None
        # start exon
        if self.cds_start is not None:
            t_pos = self.cds_start + p
        # all-coding and stop exons can be calculated the same way
        else:
            t_pos = p - self.cds_pos + self.start

        # error checking to make sure p was a proper transcript pos and inside CDS
        if self.contains_transcript_pos(t_pos) is False:
            return None
        # special case - single exon CDS
        if self.cds_stop is not None and self.cds_start is not None:
            if self.cds_start <= t_pos < self.cds_stop:
                return t_pos
        elif self.cds_stop is not None and t_pos >= self.cds_stop:
            return None
        elif self.cds_start is not None and t_pos < self.cds_start:
            return None
        else:
            return t_pos


class ChromosomeInterval(object):
    """
    Represents an interval of a chromosome. BED coordinates, strand is True,
    False or None (if no strand)
    interval arithmetic adapted from http://code.activestate.com/recipes/576816-interval/
    """
    __slots__ = ('chromosome', 'start', 'stop', 'strand')    # conserve memory

    def __init__(self, chromosome, start, stop, strand):
        self.chromosome = str(chromosome)
        self.start = int(start)    # 0 based
        self.stop = int(stop)      # exclusive
        if strand not in [True, False, None]:
            strand = convert_strand(strand)
        self.strand = strand       # True or False

    def __len__(self):
        return abs(self.stop - self.start)

    def __cmp__(self, other):
        """
        compares this chromosome interval to another
        """
        if None == other:
            return 1
        start_cmp = cmp(self.start, other.start)
        if start_cmp != 0:
            return start_cmp
        else:
            return cmp(self.stop, other.stop)

    def __contains__(self, other):
        return self.start <= other < self.stop

    def __add__(self, other):
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        return ChromosomeInterval(self.chromosome, self.start + other.start, self.stop + other.stop, self.strand)

    def __sub__(self, other):
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        return ChromosomeInterval(self.chromosome, self.start - other.start, self.stop - other.stop, self.strand)

    @property
    def is_null(self):
        if len(self) == 0:
            return True
        return False

    def intersection(self, other):
        """
        returns a new ChromosomeInterval representing the overlap between these intervals
        returns None if there is no overlap
        """
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.stop <= other.start:
            return None
        if self == other:
            return self
        if self.stop <= other.stop:
            return ChromosomeInterval(self.chromosome, other.start, self.stop, self.strand)
        else:
            return ChromosomeInterval(self.chromosome, other.start, other.stop, self.strand)

    def complement(self, size):
        """
        returns two new ChromosomeIntervals representing the complement of this interval.
        Requires a input chromosome size
        """
        return [ChromosomeInterval(self.chromosome, 0, self.start, self.strand),
                ChromosomeInterval(self.chromosome, self.stop, size, self.strand)]

    def union(self, other):
        """
        Returns either one or two ChromosomeIntervals representing the union of this interval with the other.
        If two are returned, that means the two intervals do not overlap.
        Returns None if the intervals are not on the same strand and chromosome
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self == other:
            return self
        if self.intersection(other) is not None:
            return self.hull(other)
        else:
            return [self, other]

    def hull(self, other):
        """
        Returns a new ChromosomeInterval representing the merged interval of these two, regardless of overlap
        I.E. if one interval is [1, 10) and another is [20, 30) this will return [1, 30)
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.subset(other):
            return ChromosomeInterval(self.chromosome, other.start, other.stop, self.strand)
        elif other.subset(self):
            return ChromosomeInterval(self.chromosome, self.start, self.stop, self.strand)
        return ChromosomeInterval(self.chromosome, self.start, other.stop, self.strand)

    def overlap(self, other):
        """
        Returns True if this interval overlaps other (if they are on the same strand and chromosome)
        """
        if self.chromosome != other.chromosome:
            return False
        if self > other:
            other, self = self, other
        return self.stop > other.start

    def subset(self, other):
        """
        Returns True if this interval is a subset of the other interval (if they are on the same strand and chromosome)
        """
        if self.chromosome != other.chromosome:
            return False
        return self.start >= other.start and self.stop <= other.stop

    def proper_subset(self, other):
        """
        same a subset, but only if other is entirely encased in this interval.
        """
        if self.chromosome != other.chromosome:
            return False
        return self.start > other.start and self.stop < other.stop

    def separation(self, other):
        """
        Returns a integer representing the start distance between these intervals
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.stop > other.start:
            return 0
        else:
            return other.start - self.stop

    def symmetric_separation(self, other):
        """
        Returns a pair of distances representing the symmetric distance between two intervals
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        return other.start - self.start, other.stop - self.stop

    def get_bed(self, rgb, name):
        """
        Returns BED tokens representing this interval. Requires a name and a rgb value. BED is BED12.
        """
        return [self.chromosome, self.start, self.stop, name, 0, convert_strand(self.strand), self.start, self.stop,
                rgb, 1, len(self), 0]

    def get_sequence(self, seq_dict, strand=True):
        """
        Returns the sequence for this intron in transcript orientation (reverse complement as necessary)
        If strand is False, returns the + strand regardless of transcript orientation.
        """
        if strand is False or self.strand is True:
            return seq_dict[self.chromosome][self.start:self.stop].upper()
        if self.strand is False:
            return reverse_complement(seq_dict[self.chromosome][self.start:self.stop]).upper()
        assert False

    def __repr__(self):
        return "ChromosomeInterval('{}', {}, {}, '{}')".format(self.chromosome, self.start, self.stop,
                                                               convert_strand(self.strand))


class Attribute(object):
    """
    Stores attributes from the gencode attribute file.
    """

    __slots__ = ("gene_id", "gene_name", "gene_type", "transcript_id", "transcript_type")

    def __init__(self, gene_id, gene_name, gene_type, transcript_id, transcript_type):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type


def gap_merge_intervals(intervals, gap):
    """
    Merges intervals within gap bases of each other
    """
    new_intervals = []
    for interval in intervals:
        if not new_intervals:
            new_intervals.append(copy.deepcopy(interval))
        elif interval.separation(new_intervals[-1]) <= gap:
            new_intervals[-1] = new_intervals[-1].hull(interval)
        else:
            new_intervals.append(copy.deepcopy(interval))
    return new_intervals


def interval_not_intersect_intervals(intervals, interval):
    """
    Takes a list of intervals and one other interval and determines if interval does not intersect with any of the
    intervals. returns True if this is the case.
    """
    intersections = []
    for target_interval in intervals:
        intersections.append(interval.intersection(target_interval))
    return intersections.count(None) == len(intersections)


def interval_not_within_wiggle_room_intervals(intervals, interval, wiggle_room=0):
    """
    Same thing as intervalNotIntersectIntervals but looks within wiggleRoom bases to count a valid lack of intersection
    """
    separation = [sum(interval.symmetric_separation(target_interval)) for target_interval in intervals]
    return not any([x <= 2 * wiggle_room for x in separation])  # we allow wiggle on both sides


def convert_strand(s):
    """
    Given a potential strand value, converts either from True/False/None
    to +/-/None depending on value
    """
    assert s in [True, False, None, "+", "-", "."]
    if s is True:
        return "+"
    elif s is False:
        return "-"
    elif s is None:
        return "."
    elif s == "-":
        return False
    elif s == "+":
        return True
    elif s == ".":
        return None


def complement(seq, comp=string.maketrans("ATGC", "TACG")):
    """
    given a sequence, return the complement.
    """
    return str(seq).translate(comp)


def reverse_complement(seq, comp=string.maketrans("ATGC", "TACG")):
    """
    Given a sequence, return the reverse complement.
    """
    return str(seq).translate(comp)[::-1]


_codon_table = {
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
    '': ''
    }


def codon_to_amino_acid(c):
    """
    Given a codon C, return an amino acid or ??? if codon unrecognized.
    Codons could be unrecognized due to ambiguity in IUPAC characters.
    """
    assert len(c) == 3
    if c is None:
        return None
    c = c.upper()
    if c in _codon_table:
        return _codon_table[c]
    return '?'


def find_offset(exon_frames, strand):
    """
    takes the exon_frames from a GenePredTranscript and returns the offset (how many bases of the CDS are out of frame)
    """
    frames = [x for x in exon_frames if x != -1]
    if len(frames) == 0:
        return 0
    if strand is True:
        if frames[0] == 0:
            offset = 0
        else:
            offset = 3 - frames[0]
    else:
        offset = frames[-1]
    return offset


def translate_sequence(sequence):
    """
    Translates a given DNA sequence to single-letter amino acid
    space. If the sequence is not a multiple of 3 it will be truncated
    silently.
    """
    result = []
    for i in xrange(0, len(sequence) - len(sequence) % 3, 3):
        result.append(codon_to_amino_acid(sequence[i: i + 3]))
    return "".join(result)


def read_codons(seq, offset=0, skip_last=True):
    """
    Provides an iterator that reads through a sequence one codon at a time.
    """
    l = len(seq)
    if skip_last:
        l -= 3
    for i in xrange(offset,  l - l % 3, 3):
            yield seq[i:i + 3]


def read_codons_with_position(seq, offset=0, skip_last=True):
    """
    Provides an iterator that reads through a sequence one codon at a time,
    returning both the codon and the start position in the sequence.
    """
    l = len(seq)
    if skip_last:
        l -= 3
    for i in xrange(offset, l - l % 3, 3):
            yield i, seq[i:i + 3]


def get_sequence_dict(file_path):
    """
    Returns a dictionary of fasta records.
    """
    return Fasta(file_path, as_raw=True)


def get_transcript_dict(gp_file):
    """
    Convenience function for creating transcript dictionaries
    """
    return {n: t for n, t in transcript_iterator(gp_file)}


def transcript_iterator(gp_file):
    """
    Given a path to a standard genePred file return a list of GenePredTranscript objects
    """
    with open(gp_file) as inf:
        for tokens in tokenize_stream(inf):
            t = GenePredTranscript(tokens)
            yield t.name, t


def get_transcript_attribute_dict(attribute_file):
    """
    Returns a dictionary mapping the transcript ID to an Attribute object.
    This stores all of the relevant information from the gencode attributes file.
    """
    attribute_dict = {}
    with open(attribute_file) as f:
        for line in f:
            line = line.split()
            gene_id, gene_name, gene_type, transcript_id, transcript_type = line
            attribute_dict[transcript_id] = Attribute(gene_id, gene_name, gene_type, transcript_id, transcript_type)
    return attribute_dict


def interval_to_bed(t, interval, rgb, name):
    """
    If you are turning interval objects into BED records, look here. t is a transcript object.
    Interval objects should always have start <= stop (+ strand chromosome ordering)
    """
    assert interval.stop >= interval.start, (t.name, t.chromosome)
    return [interval.chromosome, interval.start, interval.stop, name + "/" + t.name, 0, convert_strand(interval.strand),
            interval.start, interval.stop, rgb, 1, interval.stop - interval.start, 0]


def splice_intron_interval_to_bed(t, intron_interval, rgb, name):
    """
    Specific case of turning an intron interval into the first and last two bases (splice sites)
    """
    interval = intron_interval
    assert interval.stop >= interval.start, (t.name, t.chromosome)
    assert interval.stop - interval.start - 2 > 2, (t.name, t.chromosome)
    block_starts = "0,{}".format(interval.stop - interval.start - 2)
    return [interval.chromosome, interval.start, interval.stop, "/".join([name, t.name]), 0,
            convert_strand(interval.strand), interval.start, interval.stop, rgb, 2, "2,2", block_starts]


def transcript_to_bed(t, rgb, name):
    """
    Convenience function for pulling BED tokens from a Transcript object.
    """
    return t.get_bed(rgb, name)


def transcript_coordinate_to_bed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in TRANSCRIPT coordinate space and returns
    a list in BED format with the specified RGB string (128,0,0 or etc) and name.
    """
    exon_stops = [x.stop for x in t.exons]
    if t.strand is True:
        # special case - we want to slice the very last base of a exon
        # we have to do this because the last base effectively has two coordinates - the slicing coordinate
        # and the actual coordinate. This is because you slice one further than you want, I.E. x[:3] returns
        # 3 bases, but x[3] is the 4th item.
        if stop in exon_stops:
            chrom_stop = t.transcript_coordinate_to_chromosome(stop - 1) + 1
        else:
            chrom_stop = t.transcript_coordinate_to_chromosome(stop)
        chrom_start = t.transcript_coordinate_to_chromosome(start)
    else:
        if stop in exon_stops:
            chrom_start = t.transcript_coordinate_to_chromosome(stop - 1)
        else:
            chrom_start = t.transcript_coordinate_to_chromosome(stop) + 1
        chrom_stop = t.transcript_coordinate_to_chromosome(start) + 1
    assert chrom_stop >= chrom_start, (t.name, t.chromosome, start, stop, name)
    return chromosome_coordinate_to_bed(t, chrom_start, chrom_stop, rgb, name)


def cds_coordinate_to_bed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CDS coordinate space and returns
    a list in BED format with the specified RGB string (128,0,0 or etc) and name.
    """
    exon_stops = [t.transcript_coordinate_to_cds(x.stop) for x in t.exons[:-1]]
    # the last exon stop will be None because it is a slicing stop, so adjust it.
    for x in t.exons:
        if x.cds_stop is not None:
            exon_stops.append(t.transcript_coordinate_to_cds(x.cds_stop - 1) + 1)
    if t.strand is True:
        # special case - we want to slice the very last base of a exon
        # we have to do this because the last base effectively has two coordinates - the slicing coordinate
        # and the actual coordinate. This is because you slice one further than you want, I.E. x[:3] returns
        # 3 bases, but x[3] is the 4th item.
        if stop in exon_stops:
            chrom_stop = t.cds_coordinate_to_chromosome(stop - 1) + 1
        else:
            chrom_stop = t.cds_coordinate_to_chromosome(stop)
        chrom_start = t.cds_coordinate_to_chromosome(start)
    else:
        if stop in exon_stops:
            chrom_start = t.cds_coordinate_to_chromosome(stop - 1)
        else:
            chrom_start = t.cds_coordinate_to_chromosome(stop) + 1
        chrom_stop = t.cds_coordinate_to_chromosome(start) + 1
    return chromosome_coordinate_to_bed(t, chrom_start, chrom_stop, rgb, name)


def chromosome_coordinate_to_bed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CHROMOSOME coordinate space and returns
    a list in BED format with the specified RGB string and name.
    """
    assert start is not None and stop is not None, (t.name, t.chromosome, start, stop, name)
    assert stop >= start, (t.name, t.chromosome, start, stop, name)
    return t.get_bed(name=name, rgb=rgb, start_offset=start, stop_offset=stop)


def chromosome_region_to_bed(t, start, stop, rgb, name):
    """
    This is different from chromosome_coordinate_to_bed - this function will not resize the BED information
    for the input transcript, but instead be any coordinate on the chromosome.
    """
    strand = convert_strand(t.strand)
    chrom = t.chromosome
    assert start is not None and stop is not None, (t.name, start, stop, name)
    assert stop >= start, (t.name, start, stop, name)
    return [chrom, start, stop, name + "/" + t.name, 0, strand, start, stop, rgb, 1, stop - start, 0]
