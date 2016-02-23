import re

import lib.seq_lib as seq_lib
import lib.comp_ann_lib as comp_ann_lib

from src.abstract_classifier import AbstractClassifier

__author__ = "Ian Fiddes"


class StartOutOfFrame(AbstractClassifier):
    """
    Is the start out of frame, according to the ExonFrames field?
    If True, reports the first 3 bases.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        for ens_id, a in self.annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(a):
                self.classify_dict[ens_id] = 0
                continue
            # remove all -1 frames because those are UTR exons
            a_frames = [x for x in a.exon_frames if x != -1]
            if a.strand is True and a_frames[0] != 0 or a.strand is False and a_frames[-1] != 0:
                self.classify_dict[ens_id] = 1
                self.details_dict[ens_id].append(seq_lib.cds_coordinate_to_bed(a, 0, 3, self.rgb, self.column))
            else:
                self.classify_dict[ens_id] = 0
        self.dump_results_to_disk()


class BadFrame(AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least 25 codons.

    Will report a BED record of the transcript if true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        for ens_id, a in self.annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(a):
                self.classify_dict[ens_id] = 0
                continue
            if a.cds_size % 3 != 0:
                bed_rec = seq_lib.chromosome_coordinate_to_bed(a, a.thick_start, a.thick_stop, self.rgb, self.column)
                self.details_dict[ens_id].append(bed_rec)
                self.classify_dict[ens_id] = 1
            else:
                self.classify_dict[ens_id] = 0
        self.dump_results_to_disk()


class BeginStart(AbstractClassifier):
    """
    Is the first 3 bases of thick_start 'ATG'?

    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        self.get_fasta()
        for ens_id, a in self.annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(a):
                self.classify_dict[ens_id] = 0
            elif a.get_cds(self.ref_seq_dict)[:3] != "ATG":
                bed_rec = seq_lib.cds_coordinate_to_bed(a, 0, 3, self.rgb, self.column)
                self.details_dict[ens_id].append(bed_rec)
                self.classify_dict[ens_id] = 1
            else:
                self.classify_dict[ens_id] = 0
        self.dump_results_to_disk()


class EndStop(AbstractClassifier):
    """
    Are the last three bases a stop codon?
    If this is NOT true, will report a BED record of the last 3 bases.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        self.get_fasta()
        stop_codons = {'TAA', 'TGA', 'TAG'}
        for ens_id, a in self.annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(a):
                self.classify_dict[ens_id] = 0
            elif a.get_cds(self.ref_seq_dict)[-3:] not in stop_codons:
                bed_rec = seq_lib.cds_coordinate_to_bed(a, a.cds_size - 3, a.cds_size, self.rgb, self.column)
                self.details_dict[ens_id].append(bed_rec)
                self.classify_dict[ens_id] = 1
            else:
                self.classify_dict[ens_id] = 0
        self.dump_results_to_disk()


class CdsGap(AbstractClassifier):
    """
    Are any of the CDS introns too short?

    Reports a BED record for each intron interval that is too short.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, cds_filter_fn=comp_ann_lib.is_cds, mult3=False, skip_n=True):
        self.get_fasta()
        for ens_id, a in self.annotation_iterator():
            for intron in a.intron_intervals:
                is_gap = comp_ann_lib.analyze_intron_gap(a, intron, self.ref_seq_dict, cds_filter_fn, skip_n, mult3)
                if is_gap is True:
                    bed_rec = seq_lib.interval_to_bed(a, intron, self.rgb, self.column)
                    self.details_dict[ens_id].append(bed_rec)
            self.classify_dict[ens_id] = len(self.details_dict[ens_id])
        self.dump_results_to_disk()


class CdsMult3Gap(CdsGap):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, cds_filter_fn=comp_ann_lib.is_cds, mult3=True, skip_n=True):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


class UtrGap(CdsGap):
    """
    Are any UTR introns too short?

    Reports on all such introns.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, cds_filter_fn=comp_ann_lib.is_not_cds, mult3=None, skip_n=True):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


class UnknownGap(CdsGap):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, cds_filter_fn=lambda intron, t: True, mult3=None, skip_n=False):
        CdsGap.run(self, cds_filter_fn, mult3, skip_n)


class CdsNonCanonSplice(AbstractClassifier):
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

    def run(self, cds_filter_fn=comp_ann_lib.is_cds, splice_dict={"GT": "AG"}):
        self.get_fasta()
        for ens_id, a in self.annotation_iterator():
            for intron in a.intron_intervals:
                splice_is_good = comp_ann_lib.analyze_splice(intron, a, self.ref_seq_dict, cds_filter_fn, splice_dict)
                if splice_is_good is True:
                    bed_rec = seq_lib.splice_intron_interval_to_bed(a, intron, self.rgb, self.column)
                    self.details_dict[ens_id].append(bed_rec)
            self.classify_dict[ens_id] = len(self.details_dict[ens_id])
        self.dump_results_to_disk()


class CdsUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def run(self, cds_filter_fn=comp_ann_lib.is_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


class UtrNonCanonSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def run(self, cds_filter_fn=comp_ann_lib.is_not_cds, splice_dict={"GT": "AG"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


class UtrUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def run(self, cds_filter_fn=comp_ann_lib.is_not_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        CdsNonCanonSplice.run(self, cds_filter_fn, splice_dict)


class SpliceContainsUnknownBases(AbstractClassifier):
    """
    Do any of the splice junctions contain unknown bases?
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_fasta()
        for ens_id, a in self.annotation_iterator():
            for intron in a.intron_intervals:
                if comp_ann_lib.short_intron(intron) is False:
                    seq = intron.get_sequence(self.ref_seq_dict, strand=True)
                    donor, acceptor = seq[:2], seq[-2:]
                    if "N" in donor or "N" in acceptor:
                        bed_rec = seq_lib.splice_intron_interval_to_bed(a, intron, self.rgb, self.column)
                        self.details_dict[ens_id].append(bed_rec)
            self.classify_dict[ens_id] = len(self.details_dict[ens_id])
        self.dump_results_to_disk()


class InFrameStop(AbstractClassifier):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 25 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        self.get_fasta()
        for ens_id, a in self.annotation_iterator():
            cds = a.get_cds(self.ref_seq_dict)
            offset = seq_lib.find_offset(a.exon_frames, a.strand)
            for i, codon in seq_lib.read_codons_with_position(cds, offset, skip_last=True):
                amino_acid = seq_lib.codon_to_amino_acid(codon)
                if amino_acid == "*":
                    bed_rec = seq_lib.cds_coordinate_to_bed(a, i, i + 3, self.rgb, self.column)
                    self.details_dict[ens_id].append(bed_rec)
            self.classify_dict[ens_id] = len(self.details_dict[ens_id])
        self.dump_results_to_disk()


class ShortCds(AbstractClassifier):
    """
    Looks to see if this transcript has a short CDS.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        for ens_id, a in self.annotation_iterator():
            if comp_ann_lib.short_cds(a) is True and a.cds_size != 0:
                bed_rec = seq_lib.cds_coordinate_to_bed(a, 0, a.cds_size, self.rgb, self.column)
                self.details_dict[ens_id].append(bed_rec)
                self.classify_dict[ens_id] = 1
            else:
                self.classify_dict[ens_id] = 0
        self.dump_results_to_disk()


class UnknownBases(AbstractClassifier):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def make_bed_recs(self, a, s, bed_rec_fn, r=re.compile("[atgcATGC][N]+[atgcATGC]")):
        for m in re.finditer(r, s):
            yield bed_rec_fn(a, m.start() + 1, m.end() - 1, self.rgb, self.column)

    def run(self, cds=False):
        self.get_fasta()
        if cds is True:
            bed_rec_fn = seq_lib.cds_coordinate_to_bed
        else:
            bed_rec_fn = seq_lib.transcript_coordinate_to_bed
        for ens_id, a in self.annotation_iterator():
            if cds is True:
                s = a.get_cds(self.ref_seq_dict)
            else:
                s = a.get_mrna(self.ref_seq_dict)
            for bed_rec in self.make_bed_recs(a, s, bed_rec_fn):
                self.details_dict[ens_id].append(bed_rec)
            self.classify_dict[ens_id] = len(self.details_dict[ens_id])
        self.dump_results_to_disk()


class UnknownCdsBases(UnknownBases):
    def run(self, cds=True):
        UnknownBases.run(self, cds)
