"""
Classifiers in comparativeAnnotator pipeline. Broken down into 3 categories, Single Genome/Comparative/Augustus
"""
from pycbio.bio.transcripts import *
from comparativeAnnotator.lib.annotation_utils import *

__author__ = "Ian Fiddes"

colors = {'input': '219,220,222',  # grey
          'mutation': '132,35,27',  # red-ish
          'assembly': '167,206,226',  # light blue
          'alignment': '35,125,191',  # blue
          'synon': '163,116,87',  # light brown
          'nonsynon': '181,216,139',  # avocado
          'generic': '152,156,45'  # grey-yellow
          }


def StartOutOfFrame(a, ref_fasta):
    """
    Is the start out of frame, according to the ExonFrames field?
    If True, reports the first 3 bases.
    """
    # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
    if short_cds(a):
        return []
    # remove all -1 frames because those are UTR exons
    a_frames = [x for x in a.exon_frames if x != -1]
    if a.strand is True and a_frames[0] != 0 or a.strand is False and a_frames[-1] != 0:
        return cds_coordinate_to_bed(a, 0, 3, colors['alignment'], 'StartOutOfFrame')
    else:
        return []


def BadFrame(a, ref_fasta):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least 25 codons.

    Will report a BED record of the transcript if true
    """
    # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
    if short_cds(a):
        return []
    if a.cds_size % 3 != 0:
        return chromosome_coordinate_to_bed(a, a.thick_start, a.thick_stop, colors['alignment'], 'BadFrame')
    else:
        return []


def BeginStart(a, ref_fasta):
    """
    Is the first 3 bases of thick_start 'ATG'?

    Returns a BED record of the first 3 bases if this is NOT true
    """
    # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
    if short_cds(a):
        return []
    elif a.get_cds(ref_fasta)[:3] != "ATG":
        bed_rec = cds_coordinate_to_bed(a, 0, 3, rgb, column)
        details_dict[ens_id].append(bed_rec)
        classify_dict[ens_id] = 1
    else:
        return 0
dump_results_to_disk()


def EndStop(a, ref_fasta):
    """
    Are the last three bases a stop codon?
    If this is NOT true, will report a BED record of the last 3 bases.
    """
    get_fasta()


stop_codons = {'TAA', 'TGA', 'TAG'}
for ens_id, a in annotation_iterator():
    # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
    if short_cds(a):
        return 0
    elif a.get_cds(ref_seq_dict)[-3:] not in stop_codons:
        bed_rec = cds_coordinate_to_bed(a, a.cds_size - 3, a.cds_size, rgb, column)
        details_dict[ens_id].append(bed_rec)
        classify_dict[ens_id] = 1
    else:
        return 0
dump_results_to_disk()


def CdsGap(a, ref_fasta):
    """
    Are any of the CDS introns too short?

    Reports a BED record for each intron interval that is too short.
    """

    @property
    def rgb(self):
        return colors["alignment"]

    def run(self, cds_filter_fn=is_cds, mult3=False, skip_n=True):
        get_fasta()
        for ens_id, a in annotation_iterator():
            for intron in a.intron_intervals:
                is_gap = analyze_intron_gap(a, intron, ref_seq_dict, cds_filter_fn, skip_n, mult3)
                if is_gap is True:
                    bed_rec = interval_to_bed(a, intron, rgb, column)
                    details_dict[ens_id].append(bed_rec)
            classify_dict[ens_id] = len(details_dict[ens_id])
        dump_results_to_disk()


class CdsMult3Gap(CdsGap):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """

    @property
    def rgb(self):
        return colors["mutation"]

    def run(self, cds_filter_fn=is_cds, mult3=True, skip_n=True):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


class UtrGap(CdsGap):
    """
    Are any UTR introns too short?

    Reports on all such introns.
    """

    @property
    def rgb(self):
        return colors["alignment"]

    def run(self, cds_filter_fn=is_not_cds, mult3=None, skip_n=True):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


class UnknownGap(CdsGap):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.
    """

    @property
    def rgb(self):
        return colors["assembly"]

    def run(self, cds_filter_fn=lambda intron, t: True, mult3=None, skip_n=False):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


def CdsNonCanonSplice(a, ref_fasta):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    Ignores cases where the splice sites are ambiguous (contains an N)

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """

    @property
    def rgb(self):
        return colors["mutation"]

    def run(self, cds_filter_fn=is_cds, splice_dict={"GT": "AG"}):
        get_fasta()
        for ens_id, a in annotation_iterator():
            for intron in a.intron_intervals:
                splice_is_good = analyze_splice(intron, a, ref_seq_dict, cds_filter_fn, splice_dict)
                if splice_is_good is True:
                    bed_rec = splice_intron_interval_to_bed(a, intron, rgb, column)
                    details_dict[ens_id].append(bed_rec)
            classify_dict[ens_id] = len(details_dict[ens_id])
        dump_results_to_disk()


class CdsUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """

    def run(self, cds_filter_fn=is_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


class UtrNonCanonSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """

    def run(self, cds_filter_fn=is_not_cds, splice_dict={"GT": "AG"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


class UtrUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """

    def run(self, cds_filter_fn=is_not_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


def SpliceContainsUnknownBases(a, ref_fasta):
    """
    Do any of the splice junctions contain unknown bases?
    """
    get_fasta()


for ens_id, a in annotation_iterator():
    for intron in a.intron_intervals:
        if short_intron(intron) is False:
            seq = intron.get_sequence(ref_seq_dict, strand=True)
            donor, acceptor = seq[:2], seq[-2:]
            if "N" in donor or "N" in acceptor:
                bed_rec = splice_intron_interval_to_bed(a, intron, rgb, column)
                details_dict[ens_id].append(bed_rec)
    classify_dict[ens_id] = len(details_dict[ens_id])
dump_results_to_disk()


def InFrameStop(a, ref_fasta):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 25 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    get_fasta()


for ens_id, a in annotation_iterator():
    cds = a.get_cds(ref_seq_dict)
    offset = find_offset(a.exon_frames, a.strand)
    for i, codon in read_codons_with_position(cds, offset, skip_last=True):
        amino_acid = codon_to_amino_acid(codon)
        if amino_acid == "*":
            bed_rec = cds_coordinate_to_bed(a, i, i + 3, rgb, column)
            details_dict[ens_id].append(bed_rec)
    classify_dict[ens_id] = len(details_dict[ens_id])
dump_results_to_disk()


def ShortCds(a, ref_fasta):
    """
    Looks to see if this transcript has a short CDS.
    """
    for ens_id, a in annotation_iterator():
        if short_cds(a) is True and a.cds_size != 0:
            bed_rec = cds_coordinate_to_bed(a, 0, a.cds_size, rgb, column)
            details_dict[ens_id].append(bed_rec)
            classify_dict[ens_id] = 1
        else:
            return 0


dump_results_to_disk()


def UnknownBases(a, ref_fasta):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """

    @property
    def rgb(self):
        return colors["assembly"]

    def make_bed_recs(self, a, s, bed_rec_fn, r=re.compile("[atgcATGC][N]+[atgcATGC]")):
        for m in re.finditer(r, s):
            yield bed_rec_fn(a, m.start() + 1, m.end() - 1, rgb, column)

    def run(self, cds=False):
        get_fasta()
        if cds is True:
            bed_rec_fn = cds_coordinate_to_bed
        else:
            bed_rec_fn = transcript_coordinate_to_bed
        for ens_id, a in annotation_iterator():
            if cds is True:
                s = a.get_cds(ref_seq_dict)
            else:
                s = a.get_mrna(ref_seq_dict)
            for bed_rec in make_bed_recs(a, s, bed_rec_fn):
                details_dict[ens_id].append(bed_rec)
            classify_dict[ens_id] = len(details_dict[ens_id])
        dump_results_to_disk()


class UnknownCdsBases(UnknownBases):
    def run(self, cds=True):
        UnknownBases.run(self, cds)
