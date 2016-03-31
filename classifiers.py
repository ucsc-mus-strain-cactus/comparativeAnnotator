"""
Classifiers in comparativeAnnotator pipeline. Broken down into 3 categories, Single Genome/Comparative/Augustus
"""
import re
import pycbio.bio.transcripts as tx_lib
import pycbio.bio.bio as bio_lib
import comparativeAnnotator.comp_lib.annotation_utils as utils

__author__ = "Ian Fiddes"


class StartOutOfFrame(utils.AbstractClassifier):
    """
    Is the start out of frame, according to the ExonFrames field?
    If True, reports the first 3 bases.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(a):
            return []
        # remove all -1 frames because those are UTR exons
        a_frames = [x for x in a.exon_frames if x != -1]
        if a.strand is True and a_frames[0] != 0 or a.strand is False and a_frames[-1] != 0:
            return [tx_lib.cds_coordinate_to_bed(a, 0, 3, self.rgb, self.name)]
        else:
            return []


class BadFrame(utils.AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least 25 codons.
    Will report a BED record of the transcript if true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, a, ref_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(a):
            return []
        if a.cds_size % 3 != 0:
            return [tx_lib.chromosome_coordinate_to_bed(a, a.thick_start, a.thick_stop, self.rgb, self.name)]
        else:
            return []


class BeginStart(utils.AbstractClassifier):
    """
    Is the first 3 bases of thick_start 'ATG'?
    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, a, ref_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(a):
            return []
        elif a.get_cds(ref_fasta)[:3] != "ATG":
            return [tx_lib.cds_coordinate_to_bed(a, 0, 3, self.rgb, self.name)]
        else:
            return []


class EndStop(utils.AbstractClassifier):
    """
    Are the last three bases a stop codon?
    If this is NOT true, will report a BED record of the last 3 bases.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        stop_codons = {'TAA', 'TGA', 'TAG'}
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(a):
            return []
        elif a.get_cds(ref_fasta)[-3:] not in stop_codons:
            return [tx_lib.cds_coordinate_to_bed(a, a.cds_size - 3, a.cds_size, self.rgb, self.name)]
        else:
            return []


class CdsGap(utils.AbstractClassifier):
    """
    Are any of the CDS introns too short?
    Reports a BED record for each intron interval that is too short.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if len(intron) <= utils.short_intron_size:
                if len(intron) % 3 != 0 and utils.is_cds(intron, a) is True:
                    bed_rec = tx_lib.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class CdsMult3Gap(CdsGap):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if utils.short_intron(intron):
                if len(intron) % 3 == 0 and utils.is_cds(intron, a) is True:
                    bed_rec = tx_lib.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class UtrGap(CdsGap):
    """
    Are any UTR introns too short?
    Reports on all such introns.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if utils.short_intron(intron) and utils.is_cds(intron, a) is False:
                bed_rec = tx_lib.interval_to_bed(a, intron, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class UnknownGap(CdsGap):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if utils.short_intron(intron):
                if "N" in intron.get_sequence(ref_fasta):
                    bed_rec = tx_lib.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class CdsNonCanonSplice(utils.AbstractClassifier):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG
    Ignores cases where the splice sites are ambiguous (contains an N)
    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta, cds_filter_fn=utils.is_cds, splice_dict={"GT": "AG"}):
        bed_recs = []
        for intron in a.intron_intervals:
            splice_is_good = utils.analyze_splice(intron, a, ref_fasta, cds_filter_fn, splice_dict)
            if splice_is_good is True:
                bed_rec = tx_lib.splice_intron_interval_to_bed(a, intron, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class CdsUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=utils.is_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class UtrNonCanonSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=utils.is_not_cds, splice_dict={"GT": "AG"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class UtrUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=utils.is_not_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class SpliceContainsUnknownBases(utils.AbstractClassifier):
    """
    Do any of the splice junctions contain unknown bases?
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if utils.short_intron(intron) is False:
                seq = intron.get_sequence(ref_fasta, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if "N" in donor or "N" in acceptor:
                    bed_rec = tx_lib.splice_intron_interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class InFrameStop(utils.AbstractClassifier):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 25 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        cds = a.get_cds(ref_fasta)
        offset = tx_lib.find_offset(a.exon_frames, a.strand)
        for i, codon in bio_lib.read_codons_with_position(cds, offset, skip_last=True):
            amino_acid = bio_lib.codon_to_amino_acid(codon)
            if amino_acid == "*":
                bed_rec = tx_lib.cds_coordinate_to_bed(a, i, i + 3, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class ShortCds(utils.AbstractClassifier):
    """
    Looks to see if this transcript has a short CDS.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        if utils.short_cds(a) is True and a.cds_size != 0:
            bed_rec = tx_lib.cds_coordinate_to_bed(a, 0, a.cds_size, self.rgb, self.name)
            return [bed_rec]
        else:
            return []


class UnknownBases(utils.AbstractClassifier):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def make_bed_recs(self, a, s, bed_rec_fn, r=re.compile("[atgcATGC][N]+[atgcATGC]")):
        for m in re.finditer(r, s):
            yield bed_rec_fn(a, m.start() + 1, m.end() - 1, self.rgb, self.name)

    def __call__(self, a, ref_fasta, cds=False):
        if cds is True:
            bed_rec_fn = tx_lib.cds_coordinate_to_bed
            s = a.get_cds(ref_fasta)
        else:
            bed_rec_fn = tx_lib.transcript_coordinate_to_bed
            s = a.get_mrna(ref_fasta)
        bed_recs = list(self.make_bed_recs(a, s, bed_rec_fn))
        return bed_recs


class UnknownCdsBases(UnknownBases):
    def __call__(self, a, ref_fasta, cds=True):
        return UnknownBases.__call__(self, a, ref_fasta, cds)


class LongTranscript(utils.AbstractClassifier):
    """
    Is this transcript unbelievably long? Filters out poor alignments.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta, size_cutoff=1 * 10 ** 6):
        if len(a) >= size_cutoff:
            return [tx_lib.transcript_to_bed(a, self.rgb, self.name)]
        else:
            return []
