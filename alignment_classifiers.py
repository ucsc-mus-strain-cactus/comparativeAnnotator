"""
Classifiers in comparativeAnnotator pipeline. Broken down into 3 categories, Single Genome/Comparative/Augustus
"""
import luigi
from luigi.contrib import sqla
from pycbio.bio.transcripts import *
from comparativeAnnotator.lib.annotation_utils import *

__author__ = "Ian Fiddes"


class AlnExtendsOffContig(AbstractAlignmentClassifier):
    """
    Does the alignment extend off of a contig or scaffold?
    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....
    OR
    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        for aln_id, aln, t in self.alignment_transcript_iterator():
            if aln.t_start == 0 and aln.q_start != 0 or aln.t_end == aln.t_size and aln.q_end != aln.q_size:
                self.details_dict[aln_id] = seq_lib.transcript_to_bed(t, self.rgb, self.column)
                self.classify_dict[aln_id] = 1
            else:
                self.classify_dict[aln_id] = 0
        self.dump_results_to_disk()


class AlnAbutsUnknownBases(AbstractAlignmentClassifier):
    """
    Do either terminal exons abut unknown bases? (Internal exons will be classified by SpliceContainsUnknownBases)
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_fasta()
        for aln_id, t in self.transcript_iterator():
            if self.seq_dict[t.chromosome][t.start - 1] == "N":
                self.classify_dict[aln_id] = 1
                left_bed_rec = t.exon_intervals[0].get_bed(self.rgb, self.column)
                self.details_dict[aln_id].append(left_bed_rec)
            if len(t.exon_intervals) > 1 and self.seq_dict[t.chromosome][t.stop] == "N":
                self.classify_dict[aln_id] = 1
                right_bed_rec = t.exon_intervals[-1].get_bed(self.rgb, "/".join([self.column, aln_id]))
                self.details_dict[aln_id].append(right_bed_rec)
            if aln_id not in self.classify_dict:
                self.classify_dict[aln_id] = 0
        self.dump_results_to_disk()


class HasOriginalIntrons(AbstractAlignmentClassifier):
    """
    Does the alignment have all original introns? It can have more (small gaps and such), but it must have all
    original introns.

    Reports a BED for each intron that is above the minimum intron size and is not a original intron.

    Classify table reports the number of missing original introns
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, fuzz_distance=5):
        for aln_id, aln, ref_aln, t, a in self.alignment_refalignment_transcript_annotation_iterator():
            ref_starts = comp_ann_lib.fix_ref_q_starts(ref_aln)
            for intron in t.intron_intervals:
                if comp_ann_lib.is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance) is False:
                    if comp_ann_lib.short_intron(intron) is False:
                        bed_rec = seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb, self.column)
                        self.details_dict[aln_id].append(bed_rec)
            aln_starts_ends = comp_ann_lib.get_adjusted_starts_ends(t, aln)
            count = 0
            for ref_exon in a.exons[1:]:
                r = [aln_start - fuzz_distance <= ref_exon.start <= aln_end + fuzz_distance for aln_start, aln_end in
                     aln_starts_ends]
                if not any(r):
                    count += 1
            self.classify_dict[aln_id] = count
        self.dump_results_to_disk()


class CodingInsertions(AbstractAlignmentClassifier):
    """
    Does the alignment introduce insertions to the target genome?

    Reports a BED record for each such insertion

    Target insertion:

    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            for start, stop, size in comp_ann_lib.insertion_iterator(a, aln, mult3):
                if start >= t.thick_start and stop < t.thick_stop:
                    bed_rec = seq_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.column)
                    self.details_dict[aln_id].append(bed_rec)
            self.classify_dict[aln_id] = len(self.details_dict[aln_id])
        self.dump_results_to_disk()


class CodingMult3Insertions(CodingInsertions):
    """
    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingInsertions.run(self, mult3=True)


class CodingDeletions(AbstractAlignmentClassifier):
    """
    Does the alignment introduce deletions that are not a multiple of 3 to the target genome?

    Reports a BED record for each deletion. Since deletions are not visible on the target genome, just has a 0 base
    record at this position.

    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            for start, stop, size in comp_ann_lib.deletion_iterator(t, aln, mult3):
                if start >= t.thick_start and stop < t.thick_stop:
                    bed_rec = seq_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.column)
                    self.details_dict[aln_id].append(bed_rec)
            self.classify_dict[aln_id] = len(self.details_dict[aln_id])
        self.dump_results_to_disk()


class CodingMult3Deletions(CodingDeletions):
    """
    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingDeletions.run(self, mult3=True)


class FrameShift(AbstractAlignmentClassifier):
    """
    Frameshifts are caused by coding indels that are not a multiple of 3. Reports a BED entry
    spanning all blocks of coding bases that are frame-shifted.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def window_starts_stops(self, t, frame_shifts):
        indel_starts, indel_stops, spans = zip(*frame_shifts)
        # calculate cumulative frame by adding each span and taking mod 3 - zeroes imply regaining frame
        # note that this code prepends a 0 to the list, offsetting all values by 1. This is useful.
        cumulative_frame = map(lambda x: x % 3, reduce(lambda l, v: (l.append(l[-1] + v) or l), spans, [0]))
        # every start is when a zero existed in the previous spot in cumulative_frame
        windowed_starts = [x for x, y in itertools.izip(indel_starts, cumulative_frame) if
                           y == 0 or x == indel_starts[0]]
        # every stop is when a zero exists at this cumulative_frame
        windowed_stops = [x for x, y in itertools.izip(indel_stops, cumulative_frame[1:]) if y == 0]
        # sanity check
        assert any([len(windowed_starts) == len(windowed_stops), len(windowed_starts) - 1 == len(windowed_stops)])
        # now we need to fix frame and stops - if this shift extends to the end of the transcript, add that stop
        # additionally, if this is a negative strand transcript, flip starts/stops so that start is always < stop
        if len(windowed_stops) < len(windowed_starts) and t.strand is False:
            windowed_stops.append(t.thick_start)
            windowed_stops, windowed_starts = windowed_starts, windowed_stops
        elif len(windowed_stops) < len(windowed_starts):
            windowed_stops.append(t.thick_stop)
        elif t.strand is False:
            windowed_stops, windowed_starts = windowed_starts, windowed_stops
        return windowed_stops, windowed_starts

    def run(self):
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            frame_shifts = list(comp_ann_lib.frame_shift_iterator(a, t, aln))
            if len(frame_shifts) == 0:
                self.classify_dict[aln_id] = 0
                continue
            windowed_stops, windowed_starts = self.window_starts_stops(t, frame_shifts)
            for start, stop in itertools.izip(windowed_starts, windowed_stops):
                bed_rec = seq_lib.chromosome_coordinate_to_bed(t, start, stop, self.rgb, self.column)
                self.details_dict[aln_id].append(bed_rec)
            self.classify_dict[aln_id] = len(self.details_dict[aln_id])
        self.dump_results_to_disk()


class AlignmentPartialMap(AbstractAlignmentClassifier):
    """
    Does the query sequence NOT map entirely?

    a.q_size != a.q_end - a.q_start

    If so, reports the entire transcript
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        for aln_id, aln, t in self.alignment_transcript_iterator():
            if aln.q_size != aln.q_end - aln.q_start:
                bed_rec = seq_lib.transcript_to_bed(t, self.rgb, self.column)
                self.details_dict[aln_id].append(bed_rec)
                self.classify_dict[aln_id] = 1
            else:
                self.classify_dict[aln_id] = 0
        self.dump_results_to_disk()


class HasOriginalStart(AbstractAlignmentClassifier):
    """
    Does the lifted over CDS have the same 3 start bases as the original transcript?

    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            cds_positions = [t.chromosome_coordinate_to_cds(
                             aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                             for i in xrange(3)]
            if None in cds_positions:
                self.details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.column))
                self.classify_dict[aln_id] = 1
            else:
                self.classify_dict[aln_id] = 0
        self.dump_results_to_disk()


class HasOriginalStop(AbstractAlignmentClassifier):
    """
    Does the lifted over CDS have the same 3 stop bases as the original transcript?

    Returns a BED record of the last 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            cds_positions = [t.chromosome_coordinate_to_cds(
                             aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                             for i in xrange(t.cds_size - 4, t.cds_size - 1)]
            if None in cds_positions:
                self.details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.column))
                self.classify_dict[aln_id] = 1
            else:
                self.classify_dict[aln_id] = 0
        self.dump_results_to_disk()


class Nonsynonymous(AbstractAlignmentClassifier):
    """
    Do any base changes introduce nonsynonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["nonsynon"]

    def run(self, equality_test=lambda target, query: target != query):
        self.get_fasta()
        for aln_id, aln, t, a in self.alignment_transcript_annotation_iterator():
            # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
            if comp_ann_lib.short_cds(t) or comp_ann_lib.short_cds(a):
                self.classify_dict[aln_id] = 0
                continue
            for i, target_codon, query_codon in comp_ann_lib.codon_pair_iterator(a, t, aln, self.seq_dict,
                                                                                 self.ref_seq_dict):
                target_aa = seq_lib.codon_to_amino_acid(target_codon)
                query_aa = seq_lib.codon_to_amino_acid(query_codon)
                if target_codon != query_codon and equality_test(target_aa, query_aa) is True:
                    bed_rec = seq_lib.cds_coordinate_to_bed(t, i, i + 3, self.rgb, self.column)
                    self.details_dict[aln_id].append(bed_rec)
            self.classify_dict[aln_id] = len(self.details_dict[aln_id])
        self.dump_results_to_disk()


class Synonymous(Nonsynonymous):
    """
    Do any base changes introduce synonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["synon"]

    def run(self, equality_test=lambda target, query: target == query):
        Nonsynonymous.run(self, equality_test)


class Paralogy(AbstractAlignmentClassifier):
    """
    Does this transcript appear more than once in the transcript dict?
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        counts = Counter(psl_lib.remove_alignment_number(aln_id) for aln_id, aln in self.alignment_iterator())
        for aln_id, t in self.transcript_iterator():
            count = counts[psl_lib.remove_alignment_number(aln_id)] - 1
            if count > 0:
                name = self.column + "_{}_Copies".format(count)
                bed_rec = seq_lib.transcript_to_bed(t, self.rgb, name)
                self.details_dict[aln_id].append(bed_rec)
            self.classify_dict[aln_id] = count
        self.dump_results_to_disk()