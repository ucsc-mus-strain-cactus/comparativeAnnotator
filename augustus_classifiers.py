import re
from collections import Counter
import pycbio.bio.transcripts as tx_lib
import pycbio.bio.intervals as interval_lib
import comparativeAnnotator.comp_lib.annotation_utils as utils


class NotSameStrand(utils.AbstractClassifier):
    """
    Does this transcript exist on the same strand?
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t):
        if aug_t.strand != t.strand:
            return [tx_lib.transcript_to_bed(aug_t, self.rgb, self.name)]
        else:
            return []


class ExonGain(utils.AbstractClassifier):
    """
    Does the augustus version of this transcript add an exon?
    This is calculated by looking at the exon boundary intervals between the genePreds
    we expect new exons to not be within the gap-merged exons of the transMap genes.
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t):
        bed_recs = []
        aug_t_intervals = aug_t.exon_intervals
        merged_t_intervals = interval_lib.gap_merge_intervals(t.exon_intervals, gap=utils.short_intron_size)
        for interval in aug_t_intervals:
            if interval_lib.interval_not_intersect_intervals(merged_t_intervals, interval):
                bed_rec = interval.get_bed(self.rgb, "/".join([self.name, aug_t.name]))
                bed_recs.append(bed_rec)
        return bed_recs


class ExonLoss(utils.AbstractClassifier):
    """
    Does the augustus version of this transcript lose an exon?
    This is calculated by looking at the exon boundary intervals between the genePreds
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t):
        bed_recs = []
        aug_t_intervals = aug_t.exon_intervals
        merged_t_intervals = interval_lib.gap_merge_intervals(t.exon_intervals, gap=utils.short_intron_size)
        for interval in merged_t_intervals:
            if interval_lib.interval_not_intersect_intervals(aug_t_intervals, interval):
                bed_rec = interval.get_bed(self.rgb, "/".join([self.name, aug_t.name]))
                bed_recs.append(bed_rec)
        return bed_recs


class NotSimilarInternalExonBoundaries(utils.AbstractClassifier):
    """
    Does the augustus version of this transcript have the same exon boundaries, within a wiggle room range?
    Returns True if any internal exons are NOT in the same boundaries, and also reports each such splice junction
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t, short_intron_size=60, wiggle_room=50):
        merged_t_intervals = interval_lib.gap_merge_intervals(t.exon_intervals, gap=short_intron_size)
        merged_t_intervals = merged_t_intervals[1:-1]
        aug_t_intervals = aug_t.exon_intervals[1:-1]
        bed_recs = []
        for interval in merged_t_intervals:
            if interval_lib.interval_not_within_wiggle_room_intervals(aug_t_intervals, interval, wiggle_room):
                bed_rec = interval.get_bed(self.rgb, "/".join([self.name, aug_t.name]))
                bed_recs.append(bed_rec)
        return bed_recs


class NotSimilarTerminalExonBoundaries(utils.AbstractClassifier):
    """
    Does the augustus version of this transcript have the same terminal exon boundaries?
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t, short_intron_size=100, wiggle_room=200):
        merged_t_intervals = interval_lib.gap_merge_intervals(t.exon_intervals, gap=short_intron_size)
        merged_t_intervals = [merged_t_intervals[0], merged_t_intervals[-1]]
        aug_t_intervals = [aug_t.exon_intervals[0], aug_t.exon_intervals[-1]]
        bed_recs = []
        for interval in merged_t_intervals:
            if interval_lib.interval_not_within_wiggle_room_intervals(aug_t_intervals, interval, wiggle_room):
                bed_rec = interval.get_bed(self.rgb, "/".join([self.name, aug_t.name]))
                bed_recs.append(bed_rec)
        return bed_recs


class NotSameStart(utils.AbstractClassifier):
    """
    Does the augustus transcript NOT have the exact same start bases as the transMap transcript?
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t):
        bed_recs = []
        if t.thick_start != aug_t.thick_start:
            bed_rec = tx_lib.cds_coordinate_to_bed(aug_t, 0, 3, self.rgb, self.name)
            bed_recs.append(bed_rec)
        return bed_recs


class NotSameStop(utils.AbstractClassifier):
    """
    Does the augustus transcript NOT have the exact same stop bases as the transMap transcript?
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def __call__(self, aug_t, t):
        bed_recs = []
        if t.thick_stop != aug_t.thick_stop:
            s = aug_t.cds_size
            bed_rec = tx_lib.cds_coordinate_to_bed(aug_t, s - 3, s, self.rgb, self.name)
            bed_recs.append(bed_rec)
        return bed_recs


def multiple_augustus_transcripts(aug_dict):
    """
    This special non-classifier function takes the entire transcript dict and produces counts of paralogy.
    """
    r = re.compile("-[0-9]+-")
    counts = Counter("-".join(r.split(aug_t.name)) for aug_t in aug_dict)
    results = {}
    for aug_t.name, aug_t in aug_dict.iteritems():
        aug_base_id = "-".join(r.split(aug_t.name))
        count = counts[aug_base_id] - 1
        if count > 1:
            name = 'MultipleAugustusTranscripts_{}_Copies'.format(count)
            bed_rec = tx_lib.transcript_to_bed(aug_t, '128,0,0', name)
            results[aug_t.name] = [bed_rec]
        else:
            results[aug_t.name] = []
    return results
