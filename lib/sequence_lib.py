"""
Convenience library for sequence information, including 2bit and bed files.

fasta functionality was removed in favor of rapidly accessible 2bit.

Original Author: Dent Earl
Modified by Ian Fiddes
"""

import string
from itertools import izip
from math import ceil, floor

from lib.twobit import TwoBitFile, TwoBitSequence

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
    
    __slots__ = ('chromosomeInterval', 'name', 'strand', 'score', 'thickStart', 'rgb',
            'thickStop', 'start', 'stop', 'intronIntervals', 'exonIntervals', 'exons',
            'cds', 'mRna')
    
    def __init__(self, bed_tokens):
        # Text BED fields
        self.name = bed_tokens[3]
        self.strand = convertStrand(bed_tokens[5])

        # Integer BED fields
        self.score = int(bed_tokens[4])
        self.thickStart = int(bed_tokens[6])
        self.thickStop = int(bed_tokens[7])
        self.start = int(bed_tokens[1])
        self.stop = int(bed_tokens[2])
        self.rgb = map(int, bed_tokens[8].split(","))

        #interval for entire transcript including introns
        self.chromosomeInterval = ChromosomeInterval(bed_tokens[0], self.start, 
                self.stop, self.strand)

        #build chromosome intervals for exons and introns
        self.exonIntervals = self._getExonIntervals(bed_tokens)
        self.intronIntervals = self._getIntronIntervals(bed_tokens)

        #build Exons mapping transcript space coordinates to chromosome
        self.exons = self._getExons(bed_tokens)

    def __len__(self):
        return  sum(x.stop-x.start for x in self.exonIntervals)

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
            exons.append(ChromosomeInterval(chrom, start + block_start, 
                    start + block_start + block_size, strand))
        return exons

    def _getIntronIntervals(self, bed_tokens):
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
                intron = ChromosomeInterval(exon.chromosome, prevExon.stop, 
                        exon.start, exon.strand)
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
        chrom, strand = bed_tokens[0], convertStrand(bed_tokens[5])

        block_count = int(bed_tokens[9])
        block_sizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        block_starts = [int(x) for x in bed_tokens[11].split(",") if x != ""]

        ##################################################################
        ## HERE BE DRAGONS
        ## this is seriously ugly code to maintain proper mapping
        ## between coordinate spaces. See the unit tests.
        ##################################################################
        if strand is False:
            block_sizes = reversed(block_sizes)
            block_starts = reversed(block_starts)

        t_pos, cds_pos = 0, None
        for block_size, block_start in izip(block_sizes, block_starts):
            #calculate transcript relative coordinates
            this_start = t_pos
            this_stop = t_pos + block_size

            #calculate chromosome relative coordinates
            this_chrom_start = chrom_start + block_start
            this_chrom_stop = chrom_start + block_start + block_size

            #calculate transcript-relative CDS positions
            #cds_pos is pos of first coding base in CDS coordinates
            this_cds_start, this_cds_stop, this_cds_pos = None, None, None
            if strand is True:
                #special case - single exon
                if block_count == 1:
                    this_cds_pos = 0
                    this_cds_start = thick_start - this_chrom_start
                    this_cds_stop = thick_stop - this_chrom_start
                #special case - entirely non-coding
                elif thick_start == thick_stop == 0:
                    this_cds_start, this_cds_stop, this_cds_pos = None, None, None
                #is this the start codon containing exon?
                elif thick_start >= this_chrom_start and thick_start < this_chrom_stop:
                    cds_pos = this_chrom_stop - thick_start
                    this_cds_pos = 0
                    this_cds_start = this_start + thick_start - this_chrom_start
                #is this the stop codon containing exon?
                elif thick_stop > this_chrom_start and thick_stop <= this_chrom_stop:
                    this_cds_pos = cds_pos
                    cds_pos += thick_stop - this_chrom_start
                    this_cds_stop = this_stop + thick_stop - this_chrom_stop
                #is this exon all coding?
                elif this_cds_stop == None and this_cds_start == None and thick_stop >= \
                        this_chrom_stop and thick_start < this_chrom_start:
                    this_cds_pos = cds_pos
                    cds_pos += block_size
            else:
                #special case - single exon
                if block_count == 1:
                    this_cds_pos = 0
                    this_cds_start = this_chrom_stop - thick_stop
                    this_cds_stop = thick_stop - this_chrom_start + 1
                #special case - entirely non-coding
                elif thick_start == thick_stop == 0:
                    this_cds_start, this_cds_stop, this_cds_pos = None, None, None
                #is this the start codon containing exon?
                elif thick_stop > this_chrom_start and thick_stop <= this_chrom_stop:
                    cds_pos = thick_stop - this_chrom_start
                    this_cds_pos = 0
                    this_cds_start = this_start + this_chrom_stop - thick_stop
                #is this the stop codon containing exon?
                elif thick_start > this_chrom_start and thick_start < this_chrom_stop:
                    this_cds_pos = cds_pos
                    cds_pos += this_chrom_stop - thick_start
                    this_cds_stop = cds_pos + this_chrom_stop - thick_start
                #is this exon all coding?
                elif this_cds_stop == None and this_cds_start == None and thick_stop >= \
                        this_chrom_stop and thick_start < this_chrom_start:
                    this_cds_pos = cds_pos
                    cds_pos += block_size                
            exons.append( Exon(this_start, this_stop, strand, this_chrom_start, 
                    this_chrom_stop, this_cds_start, this_cds_stop, this_cds_pos) )
            t_pos += block_size

        return exons

    def getMRna(self, twoBitFileObj):
        """
        Returns the mRNA sequence for this transcript based on a TwoBitFile object.
        and the start/end positions and the exons. Sequence returned in
        5'-3' transcript orientation.
        """
        if hasattr(self, "mRna"):
            return self.mRna
        sequence = twoBitFileObj[self.chromosomeInterval.chromosome]
        assert self.chromosomeInterval.stop <= len(sequence)
        s = []
        for e in self.exonIntervals:
            s.append(sequence[e.start : e.stop])
        if self.chromosomeInterval.strand is True:
            mRna = "".join(s)
        else:
            mRna = reverseComplement("".join(s))
        self.mRna = mRna
        return mRna

    def getCds(self, twoBitFileObj):
        """
        Return the CDS sequence (as a string) for the transcript 
        (based on the exons) using a TwoBitFile object as the sequence source.
        The returned sequence is in the correct 5'-3' orientation (i.e. it has
        been reverse complemented if necessary).
        """
        if hasattr(self, "cds"):
            return self.cds
        sequence = twoBitFileObj[self.chromosomeInterval.chromosome]
        assert self.chromosomeInterval.stop <= len(sequence)
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
        if not self.chromosomeInterval.strand:
            cds = reverseComplement("".join(s))
        else:
            cds = "".join(s)
        self.cds = cds
        return cds

    def getCdsLength(self):
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

    def getProteinSequence(self, twoBitFileObj):
        """
        Returns the translated protein sequence for this transcript in single
        character space.
        """
        cds = self.getCds(twoBitFileObj)
        if len(cds) < 3:
            return ""
        return translateSequence(self.getCds(twoBitFileObj))

    def getIntronSequences(self, twoBitFileObj):
        """
        Returns a list of strings representing each intron in 5'-3' transcript
        orientation
        """
        sequence = twoBitFileObj[self.chromosomeInterval.chromosome]
        assert self.chromosomeInterval.stop <= len(sequence)
        introns = []
        prevExon = self.exonIntervals[0]
        for nextExon in self.exonIntervals[1:]:
            assert nextExon.strand == prevExon.strand
            assert nextExon.start > prevExon.stop
            start, stop = prevExon.stop, nextExon.start
            intron = sequence[start : stop]
            if self.strand is True:
                introns.append(intron)
            else:
                introns.append(reverseComplement(intron))
            prevExon = nextExon
        if self.strand is True:
            return introns
        else:
            return introns[::-1]

    def transcriptCoordinateToCds(self, p):
        """
        Takes a transcript-relative position and converts it to CDS coordinates.
        Will return None if this transcript coordinate is non-coding.
        Transcript/CDS coordinates are 0-based half open on 5'->3' transcript orientation.
        """
        for exon in self.exons:
            if exon.containsTranscriptPos(p):
                return exon.transcriptPosToCdsPos(p)

    def transcriptCoordinateToChromosome(self, p):
        """
        Takes a mRNA-relative position and converts it to chromosome position.
        Take a look at the docstring in the Exon class method chromPosToTranscriptPos
        for details on how this works.
        """
        for exon in self.exons:
            if exon.containsTranscriptPos(p):
                return exon.transcriptPosToChromPos(p)

    def chromosomeCoordinateToTranscript(self, p):
        """
        Takes a chromosome-relative position and converts it to transcript 
        coordinates. Transcript coordinates are 0-based half open on
        5'->3' transcript orientation.
        """
        for exon in self.exons:
            if exon.containsChromPos(p):
                return exon.chromPosToTranscriptPos(p)

    def chromosomeCoordinateToCds(self, p):
        """
        Takes a chromosome-relative position and converts it to CDS coordinates.
        Will return None if this chromosome coordinate is not in the CDS.
        """
        for exon in self.exons:
            if exon.containsChromPos(p):
                return exon.chromPosToCdsPos(p)

    def cdsCoordinateToTranscript(self, p):
        """
        Takes a CDS-relative position and converts it to Transcript coordinates.
        """
        for exon in self.exons:
            if exon.containsCdsPos(p):
                return exon.cdsPosToTranscriptPos(p)

    def cdsCoordinateToChromosome(self, p):
        """
        Takes a CDS-relative position and converts it to Chromosome coordinates.
        """
        for exon in self.exons:
            if exon.containsCdsPos(p):
                return exon.cdsPosToChromPos(p)

    def cdsCoordinateToAminoAcid(self, p, twoBitFileObj):
        """
        Takes a CDS-relative position and a TwoBitFile object that contains this
        transcript and returns the amino acid at that CDS position.
        Returns None if this is invalid.
        """
        cds = self.getCds(twoBitFileObj)
        if p >= len(cds) or p < 0:
            return None
        #we add 0.1 to the ceiling to make multiples of 3 work
        start, stop = int(floor(p / 3.0) * 3),  int(ceil((p + 0.1) / 3.0) * 3)
        if stop - start != 3:
            return None
        codon = cds[start : stop]
        return codonToAminoAcid(codon)

    def transcriptCoordinateToAminoAcid(self, p, twoBitFileObj):
        """
        Takes a transcript coordinate position and a TwoBitFile object that contains
        this transcript and returns the amino acid at that transcript position.
        If this position is not inside the CDS returns None.
        """
        cds_pos = self.transcriptCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, twoBitFileObj)

    def chromosomeCoordinateToAminoAcid(self, p, twoBitFileObj):
        """
        Takes a chromosome coordinate and a TwoBitFile object that contains this
        transcript and returns the amino acid at that chromosome position.
        Returns None if this position is not inside the CDS.
        """
        cds_pos = self.chromosomeCoordinateToCds(p)
        if cds_pos is None:
            return None
        return self.cdsCoordinateToAminoAcid(cds_pos, twoBitFileObj)


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
    #adding slots for cdsStartStat, cdsEndStat, exonFrames
    __slots__ = ('cdsStartStat', 'cdsEndStat', 'exonFrames')
    
    def __init__(self, gene_pred_tokens):
        # Text genePred fields
        self.name = gene_pred_tokens[0]
        self.strand = convertStrand(gene_pred_tokens[2])

        # Integer genePred fields
        self.score = 0 # no score in genePred files
        self.thickStart = int(gene_pred_tokens[5])
        self.thickStop = int(gene_pred_tokens[6])
        self.start = int(gene_pred_tokens[3])
        self.stop = int(gene_pred_tokens[4])
        self.rgb = [0, 128, 0] #no RGB in genePred files

        # genePred specific fields
        self.cdsStartStat = gene_pred_tokens[12]
        self.cdsEndStat = gene_pred_tokens[13]
        self.exonFrames = [int(x) for x in gene_pred_tokens[14].split(",") if x != ""]

        #create a fake BED entry to pass to the interval making stuff
        blockCount = gene_pred_tokens[7]
        blockStarts = [int(x) for x in gene_pred_tokens[8].split(",") if x != ""]
        blockEnds = [int(x) for x in gene_pred_tokens[9].split(",") if x != ""]
        blockSizes = [e - s for e,s in izip(blockEnds, blockStarts)]
        bed_tokens = [gene_pred_tokens[1], self.start, self.stop, self.name, self.score, gene_pred_tokens[2], 
                self.thickStart, self.thickStop, ",".join(map(str,self.rgb)), blockCount, 
                ",".join(map(str,blockSizes)), ",".join(map(str,blockStarts))]

        #interval for entire transcript including introns
        self.chromosomeInterval = ChromosomeInterval(bed_tokens[0], self.start, 
                self.stop, self.strand)

        #build chromosome intervals for exons and introns
        self.exonIntervals = self._getExonIntervals(bed_tokens)
        self.intronIntervals = self._getIntronIntervals(bed_tokens)

        #build Exons mapping transcript space coordinates to chromosome
        self.exons = self._getExons(bed_tokens)


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
    __slots__ = ('start', 'stop', 'strand', 'chromStart', 'chromStop', 
            'cdsStart', 'cdsStop', 'cdsPos')

    def __init__(self, start, stop, strand, chromStart, chromStop, cdsStart,
                cdsStop, cdsPos):
        assert chromStop - chromStart == stop - start
        self.strand = strand
        #start, stop are transcript-relative coordinates
        self.start = start
        self.stop = stop
        self.strand = strand
        self.chromStart = chromStart
        self.chromStop = chromStop
        #cdsStart/cdsStop are transcript-coordinate
        #None if not stop/start codon in this exon
        self.cdsStart = cdsStart
        self.cdsStop = cdsStop
        #cdsPos is cds coordinate
        self.cdsPos = cdsPos

    def __len__(self):
        return self.stop - self.start

    def containsChromPos(self, p):
        """does this exon contain a given chromosome position?"""
        if p is None:  return None
        elif p >= self.chromStart and p < self.chromStop:
            return True
        return False

    def containsTranscriptPos(self, p):
        """does this exon contain a given transcript position?"""
        if p is None: return None
        elif p >= self.start and p < self.stop:
            return True
        return False

    def containsCdsPos(self, p):
        """does this exon contain a given CDS position?"""
        if p is None: return None
        #see cdsPosToTranscriptPos - it will return None on invalid CDS positions
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

        (-) t.pos       4 3 2     1 0
        (+) t.pos       0 1 2     3 4
        transcript  - - A T T - - T G -
        chrom seq   G T A T T C T T G G
        chrom pos   0 1 2 3 4 5 6 7 8 9

        equivalent (negative strand) BED record:
        ['chr1', '2', '9', 'A', '0', '-', '2', '9', '0,128,0', '2', '3,2', '0,5']

        equivalent (positive strand) BED record:
        ['chr1', '2', '9', 'A', '0', '+', '2', '9', '0,128,0', '2', '3,2', '0,5']

        """
        if p is None: return None
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
        #convert to transcript space
        t_pos = self.chromPosToTranscriptPos(p)
        return self.transcriptPosToCdsPos(t_pos)

    def transcriptPosToCdsPos(self, p):
        """
        Given a transcript position, report cds-relative position.
        Returns None if this transcript position is not coding.
        Some examples:

        (+) cds.pos   0 1 2 3     4 5
        (+) t.pos   0 1 2 3 4     5 6 7 
        transcript  g T A T T - - T G g
        chrom seq   G T A T T C T T G G
        chrom pos   0 1 2 3 4 5 6 7 8 9
        equivalent BED record:
        ['chr1', '0', '10', 'A', '0', '+', '1', '9', '0,128,0', '2', '5,3', '0,7']

        (-) cds.pos   5 4 3 2     1 0
        (-) t.pos   7 6 5 4 3     2 1 0
        transcript  g T A T T - - T G g
        chrom seq   G T A T T C T T G G
        chrom pos   0 1 2 3 4 5 6 7 8 9
        equivalent BED record:
        ['chr1', '0', '10', 'A', '0', '-', '1', '9', '0,128,0', '2', '5,3', '0,7']

        (+) cds.pos   0 1   2 3     4
        (+) t.pos   0 1 2   3 4     5 6
        transcript  g T A - T C - - G g
        chrom seq   G T A T T C T T G G
        chrom pos   0 1 2 3 4 5 6 7 8 9
        equivalentt BED record:
        ['chr1', '0', '10', 'A', '0', '+', '1', '9', '0,128,0', '3', '3,2,2', '0,4,8']
        """
        if p is None: return None
        elif p < self.start or p >= self.stop:
            return None

        #is this a coding exon?
        elif self.containsCds() is False:
            return None

        #special case of single exon gene
        if self.cdsStart is not None and self.cdsStop is not None:
            if p < self.cdsStart or p >= self.cdsStop:
                return None
            return p - self.cdsStart
        #exon contains start codon
        elif self.cdsStart is not None:
            #make sure we are within CDS
            if self.cdsStart > p:
                return None
            return p - self.cdsStart
        #exon contains stop codon
        elif self.cdsStop is not None:
            #make sure we are within CDS
            if self.cdsStop <= p:
                return None
            return self.cdsPos + p - self.start
        #exon must be entirely coding
        else:
            return self.cdsPos + p - self.start

    def transcriptPosToChromPos(self, p):
        """
        Given a transcript position, returns the chromosome position
        if it is within this exon otherwise return None
        """
        if p is None: return None
        #0 based half open
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
        if p is None: return None
        #not a coding exon
        elif self.containsCds() is False:
            return None
        #start exon    
        if self.cdsStart is not None:
            t_pos = self.cdsStart + p
        #all-coding and stop exons can be calculated the same way
        else:
            t_pos = p - self.cdsPos + self.start
        
        #error checking to make sure p was a proper transcript pos and inside CDS
        if self.containsTranscriptPos(t_pos) is False:
            return None
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
        return (self.chromosome == other.chromosome and
                        self.start == other.start and
                        self.stop == other.stop and
                        self.strand == other.strand)

    def __cmp__(self, cI):
        return cmp((self.chromosome, self.start, self.stop, self.strand),
                             (cI.chromosome, cI.start, cI.stop, cI.strand))

    def __len__(self):
        return self.stop - self.start

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


class Attribute(object):
    """
    Stores attributes from the gencode attribute file.
    """
    
    __slots__ = ("geneID", "geneName", "geneType", "transcriptID", "transcriptType")
    
    def __init__(self, geneID, geneName, geneType, transcriptID, transcriptType):
        self.geneID = geneID
        self.geneName = geneName
        self.geneType = geneType
        self.transcriptID = transcriptID
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
    '': ''
    }

def codonToAminoAcid(c):
    """
    Given a codon C, return an amino acid or ??? if codon unrecognized.
    Codons could be unrecognized due to ambiguity IUPAC characters.
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
    #truncate sequence to multiple of 3
    sequence = sequence[:len(sequence) - len(sequence) % 3]
    result = []
    for i in xrange(0, len(sequence), 3):
        result.append(codonToAminoAcid(sequence[i : i + 3]))
    return "".join(result)


def readCodons(seq):
    """
    Provide an iterator that reads through a sequence one codon at a time.
    """
    for i in xrange(0, len(seq), 3):
        yield seq[i:i+3]


def readTwoBit(file_path):
    """
    Returns a dictionary that can randomly access two bit files.
    Acts as a wrapper around the TwoBitFile class in twobitreader.py.
    """
    return TwoBitFile(file_path)


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
                raise RuntimeError('transcriptListToDict: Discovered a '
                         'duplicate transcript %s %s'
                         % (t.name, t.chromosomeInterval.chromosome))
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
            geneID, geneName, geneType, geneStatus, transcriptID, transcriptName, \
                    transcriptType, transcriptStatus, havanaGeneID, havanaTranscriptID, \
                    ccdsID, level, transcriptClass = line
            attribute_dict[transcriptID] = Attribute(geneID, geneName, geneType, 
                    transcriptID, transcriptType)
    return attribute_dict