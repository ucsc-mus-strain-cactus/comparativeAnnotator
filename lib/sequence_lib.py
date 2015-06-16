"""
Convenience library for sequence information, including BED/genePred files and fasta files
Needs the python library pyfaidx installed.

Original Author: Dent Earl
Modified by Ian Fiddes
"""

import string
from itertools import izip
from math import ceil, floor
from pyfaidx import Fasta, FastaRecord, FetchError


def new_get_item(self, n):
    """
    # I have to over-write the pyfaidx way of returning sequences so that it returns upper case strings and not unicode
    """
    try:
        if isinstance(n, slice):
            start, stop, step = n.start, n.stop, n.step
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
            if stop < 0:
                stop = len(self) + stop
            if start < 0:
                start = len(self) + start
            return str(self._fa.get_seq(self.name, start + 1, stop)[::step]).upper()
        elif isinstance(n, int):
            if n < 0:
                n = len(self) + n
            return str(self._fa.get_seq(self.name, n + 1, n + 1)).upper()
    except FetchError:
        raise

FastaRecord.__getitem__ = new_get_item


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

    __slots__ = ('name', 'strand', 'score', 'thickStart', 'rgb', 'thickStop', 'start', 'stop', 'intronIntervals',
                 'exonIntervals', 'exons', 'cds', 'mRna', 'blockSizes', 'blockStarts', 'blockCount', 'chromosome',
                 'cdsSize', 'transcriptSize')

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
        self.blockSizes = bed_tokens[10]
        self.blockStarts = bed_tokens[11]
        # build chromosome intervals for exons and introns
        self.exonIntervals = self._getExonIntervals(bed_tokens)
        self.intronIntervals = self._getIntronIntervals()
        # build Exons mapping transcript space coordinates to chromosome
        self.exons = self._getExons(bed_tokens)
        # calculate sizes
        self._getCdsSize()
        self._getSize()

    def __len__(self):
        return self.transcriptSize

    def __cmp__(self, transcript):
        return cmp((self, self.name), (transcript, transcript.name))

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
            assert self.chromosomeCoordinateToTranscript(start_offset) is not None   # no intron records
            return [self.chromosome, start_offset, stop_offset, name, self.score, convertStrand(self.strand),
                    start_offset, stop_offset, rgb, 1, 0, 0]

        def _moveStart(exonIntervals, blockCount, blockStarts, blockSizes, start, start_offset):
            toRemove = len([x for x in exonIntervals if x.start <= start_offset and x.stop <= start_offset])
            assert toRemove < len(exonIntervals)
            if toRemove > 0:
                blockCount -= toRemove
                blockSizes = blockSizes[toRemove:]
                start += blockStarts[toRemove]
                new_block_starts = [0]
                for i in xrange(toRemove, len(blockStarts) - 1):
                    new_block_starts.append(blockStarts[i + 1] - blockStarts[i] + new_block_starts[-1])
                blockStarts = new_block_starts
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
            if start + blockStarts[-1] < stop_offset < stop:
                blockSizes[-1] = stop_offset - start - blockStarts[-1] 
                stop = stop_offset
            return stop, blockCount, blockStarts, blockSizes

        blockCount = int(self.blockCount)
        blockStarts = map(int, self.blockStarts.split(","))
        blockSizes = map(int, self.blockSizes.split(","))
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

    def _getExonIntervals(self, bed_tokens):
        """
        Gets a list of exon intervals in chromosome coordinate space.
        These exons are on (+) strand ordering regardless of transcript strand.
        This means (-) strand genes will be represented backwards
        """
        exons = []
        start, stop = int(bed_tokens[1]), int(bed_tokens[2])
        chrom, strand = bed_tokens[0], convertStrand(bed_tokens[5])

        block_sizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        block_starts = [int(x) for x in bed_tokens[11].split(",") if x != ""]

        for block_size, block_start in izip(block_sizes, block_starts):
            exons.append(ChromosomeInterval(chrom, start + block_start, start + block_start + block_size, strand))
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
                intron = ChromosomeInterval(exon.chromosome, prevExon.stop, exon.start, exon.strand)
                introns.append(intron)
            prevExon = exon
        return introns

    def _getExons(self, bed_tokens):
        """
        Get a list of Exons representing the exons in transcript coordinate
        space. This is in transcript order. See the Exon class for more.
        """
        exons = []
        chrom_start, chrom_stop = int(bed_tokens[1]), int(bed_tokens[2])
        thick_start, thick_stop = int(bed_tokens[6]), int(bed_tokens[7])
        if thick_start == thick_stop:
            thick_start = thick_stop = 0
        chrom, strand = bed_tokens[0], convertStrand(bed_tokens[5])

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

    def _getSize(self):
        self.transcriptSize = sum(x.stop - x.start for x in self.exonIntervals)

    def _getCdsSize(self):
        l = 0
        for e in self.exonIntervals:
            if self.thickStart < e.start and e.stop < self.thickStop:
                # squarely in the CDS
                l += e.stop - e.start
            elif e.start <= self.thickStart < e.stop < self.thickStop:
                # thickStart marks the start of the CDS
                l += e.stop - self.thickStart
            elif e.start <= self.thickStart and self.thickStop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                l += self.thickStop - self.thickStart
            elif self.thickStart < e.start < self.thickStop <= e.stop:
                # thickStop marks the end of the CDS
                l += self.thickStop - e.start
        self.cdsSize = l

    def getCdsLength(self):
        return self.cdsSize

    def getTranscriptLength(self):
        return self.transcriptSize

    def getMRna(self, seqDict):
        """
        Returns the mRNA sequence for this transcript based on a Fasta object.
        and the start/end positions and the exons. Sequence returned in
        5'-3' transcript orientation.
        """
        if hasattr(self, "mRna"):
            return self.mRna
        sequence = seqDict[self.chromosome]
        assert self.stop <= len(sequence)
        s = []
        for e in self.exonIntervals:
            s.append(sequence[e.start:e.stop])
        if self.strand is True:
            mRna = "".join(s)
        else:
            mRna = reverseComplement("".join(s))
        self.mRna = mRna
        return mRna

    def getSequence(self, seqDict):
        """
        Returns the entire chromosome sequence for this transcript, (+) strand orientation.
        """
        sequence = seqDict[self.chromosome]
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
        # make sure this isn't a non-coding gene
        if self.thickStart == self.thickStop == 0:
            return ""
        s = []
        for e in self.exonIntervals:
            if self.thickStart < e.start and e.stop < self.thickStop:
                # squarely in the CDS
                s.append(sequence[e.start:e.stop])
            elif e.start <= self.thickStart < e.stop < self.thickStop:
                # thickStart marks the start of the CDS
                s.append(sequence[self.thickStart:e.stop])
            elif e.start <= self.thickStart and self.thickStop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                s.append(sequence[self.thickStart: self.thickStop])
            elif self.thickStart < e.start < self.thickStop <= e.stop:
                # thickStop marks the end of the CDS
                s.append(sequence[e.start:self.thickStop])
        if not self.strand:
            cds = reverseComplement("".join(s)).upper()
        else:
            cds = "".join(s).upper()
        self.cds = cds
        return cds

    def getTranscriptCoordinateCdsStart(self):
        """
        Returns the transcript-relative position of the CDS start
        """
        return self.chromosomeCoordinateToTranscript(self.thickStart)

    def getTranscriptCoordinateCdsStop(self):
        """
        Returns the transcript-relative position of the CDS stop
        """
        return self.chromosomeCoordinateToTranscript(self.thickStop)

    def getChromosomeCoordinateCdsStart(self):
        """
        Returns the chromosome-relative position of the CDS start.
        This is not thickStart if a negative strand gene.
        Therefore, no gaurantee it is smaller than the stop.
        """
        if self.strand is True:
            return self.thickStart
        else:
            return self.thickStop - 1

    def getChromosomeCoordinateCdsStop(self):
        """
        Returns the chromosome-relative position of the CDS stop.
        This is not thickStart if a negative strand gene.
        Therefore, no gaurantee it is larger than the stop.
        """
        if self.strand is True:
            return self.thickStop
        else:
            return self.thickStart - 1

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
                yield intron.getSequence(seqDict)
        else:
            for intron in reversed(self.intronIntervals):
                yield intron.getSequence(seqDict)

    def getIntronSequences(self, seqDict):
        """
        Wrapper for intronSequenceIterator that returns a list of sequences.
        """
        return list(self.intronSequenceIterator(seqDict))

    def transcriptCoordinateToCds(self, p):
        """
        Takes a transcript-relative position and converts it to CDS coordinates.
        Will return None if this transcript coordinate is non-coding.
        Transcript/CDS coordinates are 0-based half open on 5'->3' transcript orientation.
        """
        for exon in self.exons:
            t = exon.transcriptPosToCdsPos(p)
            if t is not None:
                return t
        return None

    def transcriptCoordinateToChromosome(self, p):
        """
        Takes a mRNA-relative position and converts it to chromosome position.
        Take a look at the docstring in the Exon class method chromPosToTranscriptPos
        for details on how this works.
        """
        for exon in self.exons:
            t = exon.transcriptPosToChromPos(p)
            if t is not None:
                return t
        return None

    def chromosomeCoordinateToTranscript(self, p):
        """
        Takes a chromosome-relative position and converts it to transcript 
        coordinates. Transcript coordinates are 0-based half open on
        5'->3' transcript orientation.
        """
        for exon in self.exons:
            t = exon.chromPosToTranscriptPos(p)
            if t is not None:
                return t
        return None

    def chromosomeCoordinateToCds(self, p):
        """
        Takes a chromosome-relative position and converts it to CDS coordinates.
        Will return None if this chromosome coordinate is not in the CDS.
        """
        for exon in self.exons:
            t = exon.chromPosToCdsPos(p)
            if t is not None:
                return t
        return None

    def cdsCoordinateToTranscript(self, p):
        """
        Takes a CDS-relative position and converts it to Transcript coordinates.
        """
        for exon in self.exons:
            t = exon.cdsPosToTranscriptPos(p)
            if t is not None:
                return t
        return None

    def cdsCoordinateToChromosome(self, p):
        """
        Takes a CDS-relative position and converts it to Chromosome coordinates.
        """
        for exon in self.exons:
            t = exon.cdsPosToChromPos(p)
            if t is not None:
                return t
        return None

    def cdsCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a CDS-relative position and a Fasta object that contains this
        transcript and returns the amino acid at that CDS position.
        Returns None if this is invalid.
        """
        cds = self.getCds(seqDict)
        if p >= len(cds) or p < 0:
            return None
        # we add 0.1 to the ceiling to make multiples of 3 work
        start, stop = int(floor(p / 3.0) * 3),  int(ceil((p + 0.1) / 3.0) * 3)
        if stop - start != 3:
            return None
        codon = cds[start: stop]
        return codonToAminoAcid(codon)

    def transcriptCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a transcript coordinate position and a Fasta object that contains
        this transcript and returns the amino acid at that transcript position.
        If this position is not inside the CDS returns None.
        """
        cds_pos = self.transcriptCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, seqDict)

    def chromosomeCoordinateToAminoAcid(self, p, seqDict):
        """
        Takes a chromosome coordinate and a Fasta object that contains this
        transcript and returns the amino acid at that chromosome position.
        Returns None if this position is not inside the CDS.
        """
        cds_pos = self.chromosomeCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, seqDict)


class GenePredTranscript(Transcript):
    """
    Represent a transcript record from a genePred file. Stores the fields from the file
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
    # adding slots for new fields
    __slots__ = ('cdsStartStat', 'cdsEndStat', 'exonFrames')

    def __init__(self, gene_pred_tokens):
        # Text genePred fields
        self.name = gene_pred_tokens[0]
        self.chromosome = gene_pred_tokens[1]
        self.strand = convertStrand(gene_pred_tokens[2])
        # Integer genePred fields
        self.score = 0                                  # no score in genePred files
        self.thickStart = int(gene_pred_tokens[5])
        self.thickStop = int(gene_pred_tokens[6])
        self.start = int(gene_pred_tokens[3])
        self.stop = int(gene_pred_tokens[4])
        self.rgb = "128,0,0"                            # no RGB in genePred files
        # genePred specific fields
        self.cdsStartStat = gene_pred_tokens[12]
        self.cdsEndStat = gene_pred_tokens[13]
        self.exonFrames = [int(x) for x in gene_pred_tokens[14].split(",") if x != ""]
        # convert genePred format coordinates to BED-like coordinates to make intervals
        self.blockCount = gene_pred_tokens[7]
        blockStarts = [int(x) for x in gene_pred_tokens[8].split(",") if x != ""]
        blockEnds = [int(x) for x in gene_pred_tokens[9].split(",") if x != ""]
        self.blockSizes = ",".join(map(str, [e - s for e, s in izip(blockEnds, blockStarts)]))
        self.blockStarts = ",".join(map(str, [x - self.start for x in blockStarts]))
        bed_tokens = [gene_pred_tokens[1], self.start, self.stop, self.name, self.score, gene_pred_tokens[2],
                      self.thickStart, self.thickStop, self.rgb, self.blockCount,
                      self.blockSizes, self.blockStarts]
        # build chromosome intervals for exons and introns
        self.exonIntervals = self._getExonIntervals(bed_tokens)
        self.intronIntervals = self._getIntronIntervals()
        # build Exons mapping transcript space coordinates to chromosome
        self.exons = self._getExons(bed_tokens)
        # calculate sizes
        self._getCdsSize()
        self._getSize()


class Exon(object):
    """
    An Exon object stores information about one exon in both
    transcript coordinates (5'->3') and chromsome coordinates (+) strand.

    Transcript coordinates and chromosome coordinates are 0-based half open.

    Has methods to convert between the two coordinates, taking strand
    into account.

    If this exon contains the start or stop codon, contains those positions
    in transcript coordinates. If this exon is coding, contains the
    CDS-coordinate position of the first base.
    """
    __slots__ = ('start', 'stop', 'strand', 'chromStart', 'chromStop', 'cdsStart', 'cdsStop', 'cdsPos')

    def __init__(self, start, stop, strand, chromStart, chromStop, cdsStart, cdsStop, cdsPos):
        assert chromStop - chromStart == stop - start
        self.strand = strand
        # start, stop are transcript-relative coordinates
        self.start = start
        self.stop = stop
        self.strand = strand
        self.chromStart = chromStart
        self.chromStop = chromStop
        # cdsStart/cdsStop are transcript-coordinate
        # None if not stop/start codon in this exon
        self.cdsStart = cdsStart
        self.cdsStop = cdsStop
        # cdsPos is cds coordinate
        self.cdsPos = cdsPos

    def __len__(self):
        return self.stop - self.start

    def containsChromPos(self, p):
        """does this exon contain a given chromosome position?"""
        if p is None:
            return None
        elif self.chromStart <= p < self.chromStop:
            return True
        return False

    def containsTranscriptPos(self, p):
        """does this exon contain a given transcript position?"""
        if p is None:
            return None
        elif self.start <= p < self.stop:
            return True
        return False

    def containsCdsPos(self, p):
        """does this exon contain a given CDS position?"""
        if p is None:
            return None
        # see cdsPosToTranscriptPos - it will return None on invalid CDS positions
        elif self.cdsPosToTranscriptPos(p) is not None:
            return True
        return False

    def containsCds(self):
        """does this exon contain CDS?"""
        if self.cdsStart == self.cdsStop == self.cdsPos == None:
            return False
        return True

    def chromPosToTranscriptPos(self, p):
        """
        Given a chromosome position, returns the transcript position
        if it is within this exon. Otherwise, returns None
        Chromosome position is always on (+) strand
        """
        if p is None:
            return None
        elif p < self.chromStart or p >= self.chromStop:
            return None
        elif self.strand is True:
            return self.start + p - self.chromStart
        else:
            return self.start + self.chromStop - 1 - p

    def chromPosToCdsPos(self, p):
        """
        Given a chromosome position, returns the CDS position if the given
        chromosome position is in fact a CDS position on this exon.
        """
        t_pos = self.chromPosToTranscriptPos(p)
        return self.transcriptPosToCdsPos(t_pos)

    def transcriptPosToCdsPos(self, p):
        """
        Given a transcript position, report cds-relative position.
        Returns None if this transcript position is not coding.
        """
        if p is None:
            return None
        elif p < self.start or p >= self.stop:
            return None
        # is this a coding exon?
        elif self.containsCds() is False:
            return None
        # special case of single exon gene
        if self.cdsStart is not None and self.cdsStop is not None:
            if p < self.cdsStart or p >= self.cdsStop:
                return None
            return p - self.cdsStart
        # exon contains start codon
        elif self.cdsStart is not None:
            # make sure we are within CDS
            if self.cdsStart > p:
                return None
            return p - self.cdsStart
        # exon contains stop codon
        elif self.cdsStop is not None:
            # make sure we are within CDS
            if self.cdsStop <= p:
                return None
            return self.cdsPos + p - self.start
        # exon must be entirely coding
        else:
            return self.cdsPos + p - self.start

    def transcriptPosToChromPos(self, p):
        """
        Given a transcript position, returns the chromosome position
        if it is within this exon otherwise return None
        """
        if p is None:
            return None
        elif p < self.start or p >= self.stop:
            return None
        elif self.strand is True:
            return p + self.chromStart - self.start
        else:
            return self.chromStop + self.start - 1 - p

    def cdsPosToChromPos(self, p):
        """
        Given a cds position, returns the chromosome position
        if the cds position is on this exon
        """
        t_pos = self.cdsPosToTranscriptPos(p)
        return self.transcriptPosToChromPos(t_pos)

    def cdsPosToTranscriptPos(self, p):
        """
        Given a CDS position, returns the transcript position if it exists.
        Otherwise returns None.
        """
        if p is None:
            return None
        # not a coding exon
        elif self.containsCds() is False:
            return None
        # start exon
        if self.cdsStart is not None:
            t_pos = self.cdsStart + p
        # all-coding and stop exons can be calculated the same way
        else:
            t_pos = p - self.cdsPos + self.start

        # error checking to make sure p was a proper transcript pos and inside CDS
        if self.containsTranscriptPos(t_pos) is False:
            return None
        # special case - single exon CDS
        if self.cdsStop is not None and self.cdsStart is not None:
            if self.cdsStart <= t_pos < self.cdsStop:
                return t_pos
        elif self.cdsStop is not None and t_pos >= self.cdsStop:
            return None
        elif self.cdsStart is not None and t_pos < self.cdsStart:
            return None
        else:
            return t_pos


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
        return (self.chromosome == other.chromosome and self.start == other.start and self.stop == other.stop and
                self.strand == other.strand)

    def __cmp__(self, cI):
        return cmp((self.chromosome, self.start, self.stop, self.strand), (cI.chromosome, cI.start, cI.stop, cI.strand))

    def __len__(self):
        return self.stop - self.start

    def size(self):
        return self.stop - self.start

    def getBed(self, name):
        """
        Returns BED tokens representing this interval. Requires a name.
        """
        return [self.chromosome, self.start, self.stop, name, 0, convertStrand(self.strand)]

    def getSequence(self, seqDict, strand=True):
        """
        Returns the sequence for this intron in transcript orientation (reverse complement as necessary)
        If strand is False, returns the + strand regardless of transcript orientation.
        """
        if strand is False:
            return seqDict[self.chromosome][self.start:self.stop]
        if self.strand is True:
            return seqDict[self.chromosome][self.start:self.stop]
        if self.strand is False:
            return reverseComplement(seqDict[self.chromosome][self.start:self.stop])
        assert False


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
    if s is True:
        return "+"
    elif s is False:
        return "-"
    elif s is None:
        return None
    elif s == "-":
        return False
    elif s == "+":
        return True


_complement = string.maketrans("ATGC", "TACG")


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
    '': ''
    }


def codonToAminoAcid(c):
    """
    Given a codon C, return an amino acid or ??? if codon unrecognized.
    Codons could be unrecognized due to ambiguity in IUPAC characters.
    """
    if c is None:
        return None
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
    for i in xrange(0, len(sequence) - len(sequence) % 3, 3):
        result.append(codonToAminoAcid(sequence[i: i + 3]))
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


def getSequenceDict(file_path):
    """
    Returns a dictionary of fasta records.
    """
    return Fasta(file_path, as_raw=True)


def getTranscripts(bedFile):
    """
    Given a path to a standard BED return a list of Transcript objects
    """
    transcripts = []
    bedFile = open(bedFile, 'r')
    for t in transcriptIterator(bedFile):
        transcripts.append(t)
    return transcripts


def getGenePredTranscripts(gpFile):
    """
    Given a path to a standard genePred file return a list of GenePredTranscript objects
    """
    transcripts = []
    gpFile = open(gpFile, 'r')
    for t in genePredTranscriptIterator(gpFile):
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
                raise RuntimeError('transcriptListToDict: Discovered a duplicate transcript {} {}'.format(t.name,
                                                                                                          t.chromosome))
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
    for line in genePredStream:
        if line != '':
            tokens = line.rstrip().split("\t")
            yield tokens


def transcriptIterator(transcriptsBedStream):
    """
    Iterates over the transcript detailed in the bed stream producing Transcript objects.
    """
    for tokens in tokenizeBedStream(transcriptsBedStream):
        yield Transcript(tokens)


def genePredTranscriptIterator(transcriptsGpStream):
    """
    Iterates over the transcript detailed in the bed stream producing Transcript objects.
    """
    for tokens in tokenizeGenePredStream(transcriptsGpStream):
        yield GenePredTranscript(tokens)


def getTranscriptAttributeDict(attributeFile):
    """
    Returns a dictionary mapping the transcript ID to an Attribute object.
    This stores all of the relevant information from the gencode attributes file.
    """
    attribute_dict = {}
    with open(attributeFile) as f: 
        for line in f:
            line = line.split()
            geneId, geneName, geneType, transcriptId, transcriptType = line
            attribute_dict[transcriptId] = Attribute(geneId, geneName, geneType, transcriptId, transcriptType)
    return attribute_dict


def intervalToBed(t, interval, rgb, name):
    """
    If you are turning interval objects into BED records, look here. t is a transcript object.
    Interval objects should always have start <= stop (+ strand chromosome ordering)
    """
    assert interval.stop >= interval.start, (t.name, t.chromosome)
    return [interval.chromosome, interval.start, interval.stop, name + "/" + t.name, 0, convertStrand(interval.strand),
            interval.start, interval.stop, rgb, 1, interval.stop - interval.start, 0]


def spliceIntronIntervalToBed(t, intronInterval, rgb, name):
    """
    Specific case of turning an intron interval into the first and last two bases (splice sites)
    """
    interval = intronInterval
    assert interval.stop >= interval.start, (t.name, t.chromosome)
    assert interval.stop - interval.start - 2 > 2, (t.name, t.chromosome)
    blockStarts = "0,{}".format(interval.stop - interval.start - 2)
    return [interval.chromosome, interval.start, interval.stop, name + "/" + t.name, 0, convertStrand(interval.strand),
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
    exonStops = [x.stop for x in t.exons]
    if t.strand is True:
        # special case - we want to slice the very last base of a exon
        # we have to do this because the last base effectively has two coordinates - the slicing coordinate
        # and the actual coordinate. This is because you slice one further than you want, I.E. x[:3] returns
        # 3 bases, but x[3] is the 4th item.
        if stop in exonStops:
            chromStop = t.transcriptCoordinateToChromosome(stop - 1) + 1
        else:
            chromStop = t.transcriptCoordinateToChromosome(stop)
        chromStart = t.transcriptCoordinateToChromosome(start)
    else:
        if stop in exonStops:
            chromStart = t.transcriptCoordinateToChromosome(stop - 1)
        else:
            chromStart = t.transcriptCoordinateToChromosome(stop) + 1
        chromStop = t.transcriptCoordinateToChromosome(start) + 1
    assert chromStop >= chromStart, (t.name, t.chromosome, start, stop, name)
    return chromosomeCoordinateToBed(t, chromStart, chromStop, rgb, name)


def cdsCoordinateToBed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CDS coordinate space and returns
    a list in BED format with the specified RGB string (128,0,0 or etc) and name.
    """
    exonStops = [t.transcriptCoordinateToCds(x.stop) for x in t.exons[:-1]]
    # the last exonstop will be None because it is a slicing stop, so adjust it.
    for x in t.exons:
        if x.cdsStop is not None:
            exonStops.append(t.transcriptCoordinateToCds(x.cdsStop - 1) + 1)
    if t.strand is True:
        # special case - we want to slice the very last base of a exon
        # we have to do this because the last base effectively has two coordinates - the slicing coordinate
        # and the actual coordinate. This is because you slice one further than you want, I.E. x[:3] returns
        # 3 bases, but x[3] is the 4th item.
        if stop in exonStops:
            chromStop = t.cdsCoordinateToChromosome(stop - 1) + 1
        else:
            chromStop = t.cdsCoordinateToChromosome(stop)
        chromStart = t.cdsCoordinateToChromosome(start)
    else:
        if stop in exonStops:
            chromStart = t.cdsCoordinateToChromosome(stop - 1)
        else:
            chromStart = t.cdsCoordinateToChromosome(stop) + 1
        chromStop = t.cdsCoordinateToChromosome(start) + 1
    return chromosomeCoordinateToBed(t, chromStart, chromStop, rgb, name)


def chromosomeCoordinateToBed(t, start, stop, rgb, name):
    """
    Takes a transcript and start/stop coordinates in CHROMOSOME coordinate space and returns
    a list in BED format with the specified RGB string and name.
    """
    strand = convertStrand(t.strand)
    chrom = t.chromosome
    assert start is not None and stop is not None, (t.name, t.chromosome, start, stop, name)
    assert stop >= start, (t.name, t.chromosome, start, stop, name)
    return t.getBed(name=name, rgb=rgb, start_offset=start, stop_offset=stop)


def chromosomeRegionToBed(t, start, stop, rgb, name):
    """
    This is different from chromosomeCoordinateToBed - this function will not resize the BED information
    for the input transcript, but instead be any coordinate on the chromosome.
    """
    strand = convertStrand(t.strand)
    chrom = t.chromosome
    assert start is not None and stop is not None, (t.name, start, stop, name)
    assert stop >= start, (t.name, start, stop, name)
    return [chrom, start, stop, name + "/" + t.name, 0, strand, start, stop, rgb, 1, stop - start, 0]
