"""
Classifiers in comparativeAnnotator pipeline. Broken down into 3 categories, Single Genome/Comparative/Augustus
"""
import itertools
from collections import Counter
import pycbio.bio.transcripts as tx_lib
import pycbio.bio.psl as psl_lib
import pycbio.bio.bio as bio_lib
import comparativeAnnotator.comp_lib.annotation_utils as utils
from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number

__author__ = "Ian Fiddes"


class AlnExtendsOffContig(utils.AbstractClassifier):
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

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        if aln.t_start == 0 and aln.q_start != 0 or aln.t_end == aln.t_size and aln.q_end != aln.q_size:
            return [tx_lib.transcript_to_bed(t, self.rgb, self.name)]
        else:
            return []


class AlnAbutsUnknownBases(utils.AbstractClassifier):
    """
    Do either terminal exons abut unknown bases? (Internal exons will be classified by SpliceContainsUnknownBases)
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        bed_recs = []
        if tgt_fasta[t.chromosome][t.start - 1] == "N":
            left_bed_rec = t.exon_intervals[0].get_bed(self.rgb, self.name)
            bed_recs.append(left_bed_rec)
        if len(t.exon_intervals) > 1 and tgt_fasta[t.chromosome][t.stop] == "N":
            right_bed_rec = t.exon_intervals[-1].get_bed(self.rgb, "/".join([self.name, t.name]))
            bed_recs.append(right_bed_rec)
        return bed_recs


class IntroducesNewIntrons(utils.AbstractClassifier):
    """
    Reports a BED for each intron that is above the minimum intron size and is not a original intron.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, fuzz_distance=5):
        ref_starts = psl_lib.invert_q_starts(ref_aln)
        bed_recs = []
        for intron in t.intron_intervals:
            if utils.is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance) is False:
                if utils.short_intron(intron) is False:
                    bed_rec = tx_lib.splice_intron_interval_to_bed(t, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class CodingInsertions(utils.AbstractClassifier):
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

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3=False):
        bed_recs = []
        if utils.short_cds(t) or utils.short_cds(a):
            return bed_recs
        for start, stop, size in utils.insertion_iterator(a, aln, mult3):
            if start >= t.thick_start and stop < t.thick_stop:
                bed_rec = tx_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class CodingMult3Insertions(CodingInsertions):
    """
    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.
    """
    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3=True):
        return CodingInsertions.__call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3)


class CodingDeletions(utils.AbstractClassifier):
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

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3=False):
        bed_recs = []
        if utils.short_cds(t) or utils.short_cds(a):
            return bed_recs
        for start, stop, size in utils.deletion_iterator(t, aln, mult3):
            if start >= t.thick_start and stop < t.thick_stop:
                bed_rec = tx_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class CodingMult3Deletions(CodingDeletions):
    """
    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.
    """
    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3=True):
        return CodingDeletions.__call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, mult3)


class FrameShift(utils.AbstractClassifier):
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

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        bed_recs = []
        if utils.short_cds(t) or utils.short_cds(a):
            return bed_recs
        frame_shifts = list(utils.frame_shift_iterator(a, t, aln))
        if len(frame_shifts) == 0:
            return bed_recs
        windowed_stops, windowed_starts = self.window_starts_stops(t, frame_shifts)
        for start, stop in itertools.izip(windowed_starts, windowed_stops):
            bed_rec = tx_lib.chromosome_coordinate_to_bed(t, start, stop, self.rgb, self.name)
            bed_recs.append(bed_rec)
        return bed_recs


class AlignmentPartialMap(utils.AbstractClassifier):
    """
    Does the query sequence NOT map entirely?

    a.q_size != a.q_end - a.q_start

    If so, reports the entire transcript
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        if aln.q_size != aln.q_end - aln.q_start:
            return [tx_lib.transcript_to_bed(t, self.rgb, self.name)]
        return []


class HasOriginalStart(utils.AbstractClassifier):
    """
    Does the lifted over CDS have the same 3 start bases as the original transcript?

    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(t) or utils.short_cds(a):
            return []
        cds_positions = [t.chromosome_coordinate_to_cds(
                         aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                         for i in xrange(3)]
        if None in cds_positions:
            bed_rec = tx_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.name)
            return [bed_rec]
        return []


class HasOriginalStop(utils.AbstractClassifier):
    """
    Does the lifted over CDS have the same 3 stop bases as the original transcript?

    Returns a BED record of the last 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(t) or utils.short_cds(a):
            return []
        cds_positions = [t.chromosome_coordinate_to_cds(
                         aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                         for i in xrange(t.cds_size - 4, t.cds_size - 1)]
        if None in cds_positions:
            bed_rec = tx_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.name)
            return [bed_rec]
        return []


class Nonsynonymous(utils.AbstractClassifier):
    """
    Do any base changes introduce nonsynonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["nonsynon"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, equality_test=lambda target, query: target != query):
        bed_recs = []
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if utils.short_cds(t) or utils.short_cds(a):
            return bed_recs
        for i, target_codon, query_codon in utils.codon_pair_iterator(a, t, aln, tgt_fasta, ref_fasta):
            target_aa = bio_lib.codon_to_amino_acid(target_codon)
            query_aa = bio_lib.codon_to_amino_acid(query_codon)
            if target_codon != query_codon and equality_test(target_aa, query_aa) is True:
                bed_rec = tx_lib.cds_coordinate_to_bed(t, i, i + 3, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class Synonymous(Nonsynonymous):
    """
    Do any base changes introduce synonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["synon"]

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, equality_test=lambda target, query: target == query):
        return Nonsynonymous.__call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, equality_test)


def paralogy(tx_dict):
    """
    This special non-classifier function takes the entire transcript dict and produces counts of paralogy.
    """
    counts = Counter(remove_alignment_number(aln_id) for aln_id in tx_dict.iterkeys())
    results = {}
    for aln_id, t in tx_dict.iteritems():
        count = counts[remove_alignment_number(aln_id)] - 1
        if count > 0:
            name = 'Paralogy_{}_Copies'.format(count)
            bed_rec = tx_lib.transcript_to_bed(t, '128,0,0', name)
            results[aln_id] = [bed_rec]
        else:
            results[aln_id] = []
    return results
