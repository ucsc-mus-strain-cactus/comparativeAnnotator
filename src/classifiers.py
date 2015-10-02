import re
from itertools import izip
from collections import defaultdict, Counter

import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib

from src.abstract_classifier import AbstractClassifier
from lib.comp_ann_lib import deletion_iterator, insertion_iterator, frame_shift_iterator, codon_pair_iterator, \
    compare_intron_to_reference


class AlignmentAbutsLeft(AbstractClassifier):
    """
    Does the alignment ext_end off the 3' end of a scaffold?
    (regardless of transcript orientation)
    If so, reports BED of entire transcript
    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....
    Doesn't need a check for pre-existing because that doesn't matter.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_alignment_dict()
        self.get_transcript_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            if aln.strand == "+" and aln.t_start == 0 and aln.q_start != 0:
                details_dict[aln_id] = seq_lib.transcript_to_bed(self.transcript_dict[aln_id], self.rgb, self.column)
                classify_dict[aln_id] = 1
            elif aln.strand == "-" and aln.t_end == aln.t_size and aln.q_end != aln.q_size:
                details_dict[aln_id] = seq_lib.transcript_to_bed(self.transcript_dict[aln_id], self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class AlignmentAbutsRight(AbstractClassifier):
    """
    Does the alignment ext_end off the 3' end of a scaffold?
    (regardless of transcript orientation)
    If so, reports BED of entire transcript
    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|
    Doesn't need a check for pre-existing because that doesn't matter.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_alignment_dict()
        self.get_transcript_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            if aln.strand == "+" and aln.t_end == aln.t_size and aln.q_end != aln.q_size:
                details_dict[aln_id] = seq_lib.transcript_to_bed(self.transcript_dict[aln_id], self.rgb, self.column)
                classify_dict[aln_id] = 1
            elif aln.strand == "-" and aln.t_start == 0 and aln.q_start != 0:
                details_dict[aln_id] = seq_lib.transcript_to_bed(self.transcript_dict[aln_id], self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class AlignmentAbutsUnknownBases(AbstractClassifier):
    """
    Do any of the exons in this alignment abut unknown bases? Defined as there being a N within 2bp of a exon
    Ignores introns below short_intron_size (which will be picked up by UnknownGap)
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, distance=2, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            intervals = [[t.exon_intervals[0].start - distance, t.exon_intervals[0].start]]
            for intron in t.intron_intervals:
                if len(intron) > short_intron_size:
                    intervals.append([intron.start, intron.start + distance])
            intervals.append([t.exon_intervals[-1].stop, t.exon_intervals[-1].stop + distance])
            for start, stop in intervals:
                seq = self.fasta[t.chromosome][start:stop].upper()
                if "N" in seq:
                    classify_dict[aln_id] = 1
                    details_dict[aln_id].append(t.get_bed(self.rgb, self.column))
                    break
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class HasOriginalIntrons(AbstractClassifier):
    """
    Does the alignment have all original introns? It can have more (small gaps and such), but it must have all
    original introns.

    Reports a BED for each intron that is above the minimum intron size and is not a original intron.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, short_intron_size=30):
        self.get_annotation_dict()
        self.get_transcript_dict()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            original_introns = {(x.start, x.stop) for x in a.intron_intervals}
            target_introns = set()
            target_intron_mapping = {}
            for intron in t.intron_intervals:
                a_start = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.start - 1)) + 1
                a_stop = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.stop))
                target_introns.add((a_start, a_stop))
                target_intron_mapping[(a_start, a_stop)] = intron
            missing_introns = original_introns - target_introns
            if len(missing_introns) != 0:
                classify_dict[aln_id] = 1
                not_original_introns = target_introns - original_introns
                for a_start, a_stop in not_original_introns:
                    intron = target_intron_mapping[(a_start, a_stop)]
                    if len(intron) >= short_intron_size:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb,
                                                                                          self.column))
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CodingInsertions(AbstractClassifier):
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
        self.get_annotation_dict()
        self.get_alignment_dict()
        self.get_transcript_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            insertions = [seq_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.column) for start, stop, size
                          in insertion_iterator(a, t, aln, mult3) if start >= t.thick_start and stop < t.thick_stop]
            if len(insertions) > 0:
                details_dict[aln_id] = insertions
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CodingMult3Insertions(CodingInsertions):
    """
    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingInsertions.run(self, mult3=True)


class CodingDeletions(AbstractClassifier):
    """
    Does the alignment introduce deletions that are not a multiple of 3 to the target genome?

    Reports a BED record for each deletion. Since deletions are not visible on the target genome, just has a 0 base
    record at this position.

    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, mult3=False):
        self.get_alignment_dict()
        self.get_transcript_dict()
        self.get_annotation_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            deletions = [seq_lib.chromosome_region_to_bed(t, start, stop, self.rgb, self.column) for start, stop, size
                         in deletion_iterator(a, t, aln, mult3) if start >= t.thick_start and stop < t.thick_stop]
            if len(deletions) > 0:
                details_dict[aln_id] = deletions
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CodingMult3Deletions(CodingDeletions):
    """
    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.
    """
    def run(self):
        CodingDeletions.run(self, mult3=True)


class StartOutOfFrame(AbstractClassifier):
    """
    StartOutOfFrame are caused when the starting CDS base of the lifted over transcript is not in the original frame.
    If True, reports the first 3 bases.

    TODO: should be flagged as previously existing in reference
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        self.get_transcript_dict()
        self.get_annotation_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codon
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            # is this is a problem in the reference?
            # remove all -1 frames because those are UTR exons
            a_frames = [x for x in a.exon_frames if x != -1]
            if a.strand is True and a_frames[0] != 0 or a.strand is False and a_frames[-1] != 0:
                classify_dict[aln_id] = 1
                details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, self.colors["input"], self.column)
                continue
            # remove all -1 frames because those are UTR exons
            t_frames = [x for x in t.exon_frames if x != -1]
            if t.strand is True and t_frames[0] != 0 or t.strand is False and t_frames[-1] != 0:
                classify_dict[aln_id] = 1
                details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.column)
                continue
            classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class FrameShift(AbstractClassifier):
    """
    Frameshifts are caused by coding indels that are not a multiple of 3. Reports a BED entry
    spanning all blocks of coding bases that are frame-shifted. Must have at least 25 codons.

    Doesn't need a check for existing in the reference because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        self.get_alignment_dict()
        self.get_transcript_dict()
        self.get_annotation_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            # do not include noncoding transcripts or lift-overs that contain less than 1 codon
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            frame_shifts = list(frame_shift_iterator(a, t, aln))
            if len(frame_shifts) == 0:
                classify_dict[aln_id] = 0
                continue
            indel_starts, indel_stops, spans = zip(*frame_shifts)
            # calculate cumulative frame by adding each span and taking mod 3 - zeroes imply regaining frame
            # note that this code prepends a 0 to the list, offsetting all values by 1. This is useful.
            cumulative_frame = map(lambda x: x % 3, reduce(lambda l, v: (l.append(l[-1] + v) or l), spans, [0]))
            # every start is when a zero existed in the previous spot in cumulative_frame
            windowed_starts = [x for x, y in izip(indel_starts, cumulative_frame) if y == 0 or x == indel_starts[0]]
            # every stop is when a zero exists at this cumulative_frame
            windowed_stops = [x for x, y in izip(indel_stops, cumulative_frame[1:]) if y == 0]
            # sanity check
            assert any([len(windowed_starts) == len(windowed_stops), len(windowed_starts) - 1 == len(windowed_stops)]),\
                (self.genome, self.column, aln_id)
            # now we need to fix frame and stops - if this shift extends to the end of the transcript, add that stop
            # additionally, if this is a negative strand transcript, flip starts/stops so that start is always < stop
            if len(windowed_stops) < len(windowed_starts) and t.strand is False:
                windowed_stops.append(t.thick_start)
                windowed_stops, windowed_starts = windowed_starts, windowed_stops
            elif len(windowed_stops) < len(windowed_starts):
                windowed_stops.append(t.thick_stop)
            elif t.strand is False:
                windowed_stops, windowed_starts = windowed_starts, windowed_stops
            for start, stop in izip(windowed_starts, windowed_stops):
                details_dict[aln_id].append(seq_lib.chromosome_coordinate_to_bed(t, start, stop, self.rgb, self.column))
            classify_dict[aln_id] = 1
        self.dump_results_to_disk(classify_dict, details_dict)


class AlignmentPartialMap(AbstractClassifier):
    """
    Does the query sequence NOT map entirely?

    a.q_size != a.q_end - a.qStart

    If so, reports the entire transcript

    Doesn't need a check for pre-existing because that is impossible.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_alignment_dict()
        self.get_transcript_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            if aln.q_size != aln.q_end - aln.q_start:
                details_dict[aln_id] = seq_lib.transcript_to_bed(self.transcript_dict[aln_id], self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class BadFrame(AbstractClassifier):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least 25 codons.

    Will report a BED record of the transcript if true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        self.get_alignment_dict()
        self.get_transcript_dict()
        self.get_annotation_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            if t.get_cds_length() % 3 != 0 and a.get_cds_length() % 3 != 0:
                details_dict[aln_id] = seq_lib.chromosome_coordinate_to_bed(t, t.thick_start, t.thick_stop, 
                                                                            self.colors["input"], self.column)
                classify_dict[aln_id] = 1
            elif t.get_cds_length() % 3 != 0:
                details_dict[aln_id] = seq_lib.chromosome_coordinate_to_bed(t, t.thick_start, t.thick_stop, self.rgb,
                                                                            self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class BeginStart(AbstractClassifier):
    """
    Does the lifted over CDS have the same 3 start bases as the original transcript?
    AND are these bases 'ATG'? Must have at least 25 codons.

    Returns a BED record of the first 3 bases if this is NOT true
    """
    @property
    def rgb(self):
        return self.colors["generic"]

    def run(self):
        self.get_transcript_dict()
        self.get_alignment_dict()
        self.get_annotation_dict()
        self.get_fasta()
        self.get_ref_fasta()
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            # do not include noncoding transcripts or lift-overs that contain less than 25 codons
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            cds_positions = [t.chromosome_coordinate_to_cds(
                             aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                             for i in xrange(3)]
            if None in cds_positions:
                details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.column)
                classify_dict[aln_id] = 1
            elif t.get_cds(self.fasta)[:3] != "ATG":
                if a.get_cds(self.ref_fasta)[:3] != "ATG":
                    details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, self.colors["input"], self.column)
                else:
                    details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CdsGap(AbstractClassifier):
    """
    Are any of the CDS introns too short? Too short default is 30 bases.

    Reports a BED record for each intron interval that is too short.

    Not currently implementing checking for pre-existing because its so rare (none in VM4)
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            for i, intron in enumerate(t.intron_intervals):
                if len(intron) >= short_intron_size:
                    continue
                elif "N" in intron.get_sequence(self.fasta):
                    continue
                elif not (intron.start >= t.thick_start and intron.stop < t.thick_stop):
                    continue
                details_dict[aln_id].append(seq_lib.interval_to_bed(t, intron, self.rgb, self.column))
                classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CdsMult3Gap(AbstractClassifier):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            for i, intron in enumerate(t.intron_intervals):
                if len(intron) >= short_intron_size:
                    continue
                elif len(intron) % 3 != 0:
                    continue
                elif "N" in intron.get_sequence(self.fasta):
                    continue
                elif not (intron.start >= t.thick_start and intron.stop < t.thick_stop):
                    continue
                details_dict[aln_id].append(seq_lib.interval_to_bed(t, intron, self.rgb, self.column))
                classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class UtrGap(AbstractClassifier):
    """
    Are any UTR introns too short? Too short is defined as less than 30bp

    Reports on all such introns.

    Not currently implementing checking for pre-existing because its so rare (397 in VM4)
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            for i, intron in enumerate(t.intron_intervals):
                if len(intron) >= short_intron_size:
                    continue
                elif "N" in intron.get_sequence(self.fasta):
                    continue
                elif intron.start >= t.thick_start and intron.stop < t.thick_stop:
                    continue
                details_dict[aln_id].append(seq_lib.interval_to_bed(t, intron, self.rgb, self.column))
                classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class UnknownGap(AbstractClassifier):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.

    Not implementing looking for pre-existing because that doesn't make sense.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            for i, intron in enumerate(t.intron_intervals):
                if len(intron) >= short_intron_size:
                    continue
                elif "N" not in intron.get_sequence(self.fasta):
                    continue
                details_dict[aln_id].append(seq_lib.interval_to_bed(t, intron, self.rgb, self.column))
                classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CdsNonCanonSplice(AbstractClassifier):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    Reports two BED records of the four offending bases.

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    canonical = {"GT": "AG"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_annotation_dict()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            for intron in t.intron_intervals:
                if len(intron) <= short_intron_size:
                    continue
                elif not (intron.start >= t.thick_start and intron.stop < t.thick_stop):
                    continue
                seq = intron.get_sequence(self.fasta, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.canonical or self.canonical[donor] != acceptor:
                    classify_dict[aln_id] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compare_intron_to_reference(intron, a, aln, self.canonical, self.ref_fasta) is True:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron,
                                                                                          self.colors["input"],
                                                                                          self.column))
                    else:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb,
                                                                                          self.column))
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class CdsUnknownSplice(AbstractClassifier):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_annotation_dict()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            for intron in t.intron_intervals:
                if len(intron) <= short_intron_size:
                    continue
                elif not (intron.start >= t.thick_start and intron.stop < t.thick_stop):
                    continue
                seq = intron.get_sequence(self.fasta, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.non_canonical or self.non_canonical[donor] != acceptor:
                    classify_dict[aln_id] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compare_intron_to_reference(intron, a, aln, self.non_canonical, self.ref_fasta) is True:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron,
                                                                                          self.colors["input"],
                                                                                          self.column))
                    else:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb, 
                                                                                          self.column))
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class UtrNonCanonSplice(AbstractClassifier):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    canonical = {"GT": "AG"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_annotation_dict()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            for intron in t.intron_intervals:
                if len(intron) <= short_intron_size:
                    continue
                elif intron.start >= t.thick_start and intron.stop < t.thick_stop:
                    continue
                seq = intron.get_sequence(self.fasta, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.canonical or self.canonical[donor] != acceptor:
                    classify_dict[aln_id] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compare_intron_to_reference(intron, a, aln, self.canonical, self.ref_fasta) is True:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, 
                                                                                          self.colors["input"],
                                                                                          self.column))
                    else:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb, 
                                                                                          self.column))
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class UtrUnknownSplice(AbstractClassifier):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}

    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self, short_intron_size=30):
        self.get_transcript_dict()
        self.get_fasta()
        self.get_alignment_dict()
        self.get_annotation_dict()
        self.get_ref_fasta()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            for intron in t.intron_intervals:
                if len(intron) <= short_intron_size:
                    continue
                elif intron.start >= t.thick_start and intron.stop < t.thick_stop:
                    continue
                seq = intron.get_sequence(self.fasta, strand=True)
                donor, acceptor = seq[:2], seq[-2:]
                if donor not in self.non_canonical or self.non_canonical[donor] != acceptor:
                    classify_dict[aln_id] = 1
                    # is this a intron that exists in the reference that also has this problem?
                    if compare_intron_to_reference(intron, a, aln, self.non_canonical, self.ref_fasta) is True:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron,
                                                                                          self.colors["input"],
                                                                                          self.column))
                    else:
                        details_dict[aln_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, self.rgb, 
                                                                                          self.column))
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class EndStop(AbstractClassifier):
    """
    Looks at the last three bases of the annotated transcript and determines if they properly match the last 3 bases
    of the lifted over transcript AND that those bases are in ('TAA', 'TGA', 'TAG'). Must have at least 25 codons.

    If this is NOT true, will report a BED record of the last 3 bases.

    Value will be NULL if there is insufficient information, which is defined as:
        1) thick_stop - thick_start <= 9: (no useful CDS annotation)
        2) this alignment was not trans-mapped
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self):
        stop_codons = ('TAA', 'TGA', 'TAG')
        self.get_alignment_dict()
        self.get_transcript_dict()
        self.get_annotation_dict()
        self.get_fasta()
        self.get_ref_fasta()
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            s = t.get_cds_length()
            cds_positions = [t.chromosome_coordinate_to_cds(
                             aln.query_coordinate_to_target(a.cds_coordinate_to_transcript(i)))
                             for i in xrange(s - 4, s - 1)]
            if None in cds_positions or t.get_cds(self.fasta)[-3:] not in stop_codons:
                # does this problem exist in the reference?
                if a.get_cds(self.ref_fasta)[-3:] not in stop_codons:
                    details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, s - 3, s, self.colors["input"], self.column)
                else:
                    details_dict[aln_id] = seq_lib.cds_coordinate_to_bed(t, s - 3, s, self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


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
        self.get_transcript_dict()
        self.get_annotation_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = self.alignment_dict[aln_id]
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            # TODO: this will miss an inframe stop if it is the last 3 bases that are not the annotated stop.
            # use the logic from EndStop to flag this
            codons = list(codon_pair_iterator(a, t, aln, self.fasta, self.ref_fasta))[:-1]
            for i, target_codon, query_codon in codons:
                if seq_lib.codon_to_amino_acid(target_codon) == "*":
                    if target_codon == query_codon:
                        details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, i, i + 3, self.colors["input"],
                                                                                  self.column))
                    else:
                        details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, i, i + 3, self.rgb, self.column))
                    classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class ShortCds(AbstractClassifier):
    """
    Looks to see if this transcript actually has a CDS, which is defined as having a CDS region of
    at least 25 codons. Adjusting cds_cutoff can change this.

    If True, reports entire transcript.

    Only reports if the original transcript had a CDS.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def run(self, cds_cutoff=75):
        self.get_transcript_dict()
        self.get_annotation_dict()
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            # do not include noncoding transcripts
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            if a.get_cds_length() < 3:
                continue
            elif a.get_cds_length() <= cds_cutoff:
                details_dict[aln_id] = seq_lib.transcript_to_bed(t, self.colors["input"], self.column)
                classify_dict[aln_id] = 1
            elif t.get_cds_length() <= cds_cutoff:
                details_dict[aln_id] = seq_lib.transcript_to_bed(t, self.rgb, self.column)
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class ScaffoldGap(AbstractClassifier):
    """
    Does this alignment span a scaffold gap? (Defined as a run of Ns more than 10bp)

    Only true if the scaffold gap is within the alignment
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self):
        self.get_fasta()
        self.get_transcript_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        r = re.compile("[ATGC][N]{11,}[ATGC]")
        for aln_id, t in self.transcript_dict.iteritems():
            for exon in t.exon_intervals:
                exon_seq = exon.get_sequence(self.fasta, strand=False)
                if r.match(exon_seq):
                    classify_dict[aln_id] = 1
                    details_dict[aln_id].append(exon.get_bed(self.rgb, self.column))
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class UnknownBases(AbstractClassifier):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def run(self, cds=False):
        self.get_transcript_dict()
        self.get_fasta()
        details_dict = {}
        classify_dict = {}
        r = re.compile("[atgcATGC][N]+[atgcATGC]")
        for aln_id, t in self.transcript_dict.iteritems():
            if cds is True:
                s = t.get_cds(self.fasta)
                tmp = [seq_lib.cds_coordinate_to_bed(t, m.start() + 1, m.end() - 1, self.rgb, self.column) for m in
                       re.finditer(r, s)]
            else:
                s = t.get_mrna(self.fasta)
                tmp = [seq_lib.transcript_coordinate_to_bed(t, m.start() + 1, m.end() - 1, self.rgb, self.column) for m
                       in re.finditer(r, s)]
            if len(tmp) > 0:
                details_dict[aln_id] = tmp
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class Nonsynonymous(AbstractClassifier):
    """
    Do any base changes introduce nonsynonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["nonsynon"]

    def run(self):
        self.get_transcript_dict()
        self.get_annotation_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            for i, target_codon, query_codon in codon_pair_iterator(a, t, aln, self.fasta, self.ref_fasta):
                if ("N" not in target_codon and target_codon != query_codon and
                        seq_lib.codon_to_amino_acid(target_codon) != seq_lib.codon_to_amino_acid(query_codon)):
                    details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, i, i + 3, self.rgb, self.column))
                    classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class Synonymous(AbstractClassifier):
    """
    Do any base changes introduce synonymous changes? Only looks at aligned pairs of codons in the frame
    of the reference annotation.
    """
    @property
    def rgb(self):
        return self.colors["synon"]

    def run(self):
        self.get_transcript_dict()
        self.get_annotation_dict()
        self.get_fasta()
        self.get_ref_fasta()
        self.get_alignment_dict()
        details_dict = defaultdict(list)
        classify_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            if aln_id not in self.transcript_dict:
                continue
            t = self.transcript_dict[aln_id]
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
                continue
            for i, target_codon, query_codon in codon_pair_iterator(a, t, aln, self.fasta, self.ref_fasta):
                if (target_codon != query_codon and
                            seq_lib.codon_to_amino_acid(target_codon) == seq_lib.codon_to_amino_acid(query_codon)):
                    details_dict[aln_id].append(seq_lib.cds_coordinate_to_bed(t, i, i + 3, self.rgb, self.column))
                    classify_dict[aln_id] = 1
            if aln_id not in classify_dict:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)


class Paralogy(AbstractClassifier):
    """
    Does this transcript appear more than once in the transcript dict?
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def run(self):
        self.get_transcript_dict()
        counts = Counter(psl_lib.remove_alignment_number(aln_id) for aln_id in self.transcript_dict)
        details_dict = {}
        classify_dict = {}
        for aln_id, t in self.transcript_dict.iteritems():
            if counts[psl_lib.remove_alignment_number(aln_id)] > 1:
                n = self.column + "_{}_Copies".format(counts[psl_lib.remove_alignment_number(aln_id)] - 1)
                b = seq_lib.transcript_to_bed(t, self.rgb, n)
                details_dict[aln_id] = b
                classify_dict[aln_id] = 1
            else:
                classify_dict[aln_id] = 0
        self.dump_results_to_disk(classify_dict, details_dict)
