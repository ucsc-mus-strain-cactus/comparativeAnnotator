"""
This file contains helper functions for comparativeAnnotator.
"""
import lib.seq_lib as seq_lib


__author__ = "Ian Fiddes"

short_intron_size = 30
short_cds_size = 75


def short_cds(t):
    """
    Many classifiers do not apply to CDS below a cutoff size.
    """
    return True if t.cds_size <= short_cds_size else False


def short_intron(intron):
    """
    Many classifiers rely on analyzing introns either above or below a cutoff size
    """
    return True if len(intron) <= short_intron_size else False


def is_cds(intron, t):
    return not (intron.start >= t.thick_start and intron.stop < t.thick_stop)


def is_not_cds(intron, t):
    return intron.start >= t.thick_start and intron.stop < t.thick_stop


def analyze_intron_gap(t, intron, seq_dict, cds_fn, skip_n, mult3):
    if skip_n is True and "N" in intron.get_sequence(seq_dict):
        return False
    elif skip_n is False and "N" not in intron.get_sequence(seq_dict):
        return False
    elif cds_fn(intron, t) is True:
        return False
    elif mult3 is True and len(intron) % 3 != 0:
        return False
    elif mult3 is False and len(intron) % 3 == 0:
        return False
    elif short_intron(intron):
        return True
    else:
        return False


def analyze_splice(intron, t, seq_dict, cds_fn, splice_sites):
    if short_intron(intron) is True:
        return False
    seq = intron.get_sequence(seq_dict, strand=True)
    donor, acceptor = seq[:2], seq[-2:]
    if cds_fn(intron, t) is True:
        return False
    elif "N" in donor or "N" in acceptor:
        return False
    elif donor not in splice_sites or splice_sites[donor] != acceptor:
        return True


def insertion_iterator(a, aln, mult3=None):
    """
    Target insertion:
    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA

    Analyze a given annotation transcript and alignment for target insertions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported. Set to None to report all.
    """
    prev_target_i = None
    exon_starts = [x.start for x in a.exons]
    for query_i in xrange(len(a)):
        if query_i in exon_starts:
            # don't call a intron an insertion
            prev_target_i = None
            continue
        target_i = aln.query_coordinate_to_target(query_i)
        if target_i is None:
            # deletion; ignore
            prev_target_i = target_i
            continue
        if prev_target_i is not None and abs(target_i - prev_target_i) != 1:
            # jumped over a insertion
            insert_size = abs(target_i - prev_target_i) - 1
            start = min(prev_target_i, target_i) + 1
            stop = max(prev_target_i, target_i)
            if mult3 is True and insert_size % 3 == 0:
                yield start, stop, insert_size
            elif mult3 is False and insert_size % 3 != 0:
                yield start, stop, insert_size
            elif mult3 is None:
                yield start, stop, insert_size
        prev_target_i = target_i


def deletion_iterator(t, aln, mult3=None):
    """
    Target deletion:
    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA

    Analyze a given transcript and alignment for target deletions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported. Set to None to report all.
    """
    prev_query_i = None
    for target_i in xrange(len(t)):
        target_chrom_i = t.transcript_coordinate_to_chromosome(target_i)
        query_i = aln.target_coordinate_to_query(target_chrom_i)
        if query_i is None:
            # insertion; ignore
            prev_query_i = query_i
            continue
        if prev_query_i is not None and abs(query_i - prev_query_i) != 1:
            # jumped over a deletion
            delete_size = abs(query_i - prev_query_i) - 1
            if t.strand is True:
                start = stop = target_chrom_i - 1
            else:
                start = stop = target_chrom_i + 1
            if mult3 is True and delete_size % 3 == 0:
                yield start, stop, -delete_size
            elif mult3 is False and delete_size % 3 != 0:
                yield start, stop, -delete_size
            elif mult3 is None:
                yield start, stop, -delete_size
        prev_query_i = query_i


def start_out_of_frame(t):
    """
    Determines if a start is out of frame by making use of the exonFrames field in a genePred.
    """
    # can't use this classifier on BED records due to missing information
    assert isinstance(t, seq_lib.GenePredTranscript)
    t_frames = [x for x in t.exon_frames if x != -1]
    return True if t.strand is True and t_frames[0] != 0 or t.strand is False and t_frames[-1] != 0 else False


def frame_shift_iterator(a, t, aln):
    """
    Yields frameshift-causing mutations. These are defined as non mult3 indels within CDS.
    """
    deletions = list(deletion_iterator(t, aln, mult3=False))
    insertions = list(insertion_iterator(a, aln, mult3=False))
    # need to fix case where indel immediately precedes indel by merging overlapping ranges and combining spans
    merged = []
    for higher in sorted(deletions + insertions):
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound, merged[-1][2] + higher[2])
            else:
                merged.append(higher)
    if t.strand is False:
        merged = [[y, x, z] for x, y, z in reversed(merged) if z % 3 != 0]
    else:
        merged = [[x, y, z] for x, y, z in merged if z % 3 != 0]
    for start, stop, span in merged:
        if start >= t.thick_start and stop < t.thick_stop:
            yield start, stop, span


def codon_pair_iterator(a, t, aln, target_seq_dict, query_seq_dict):
    """
    Inputs:
    Transcript objects representing the annotation (query) transcript and the target transcript.
    PslRow object that represents the alignment between the transcript objects.
    pyfaidx Fasta objects that contain the genomic sequence for these two transcripts

    Order is (target_cds_pos, target, query)
    """
    target_cds = t.get_cds(target_seq_dict)
    query_cds = a.get_cds(query_seq_dict)
    a_frames = [x for x in a.exon_frames if x != -1]
    a_offset = seq_lib.find_offset(a_frames, a.strand)
    for i in xrange(a_offset, a.cds_size - a.cds_size % 3, 3):
        target_cds_positions = [t.chromosome_coordinate_to_cds(
                                aln.query_coordinate_to_target(
                                a.cds_coordinate_to_transcript(j)))
                                for j in xrange(i, i + 3)]
        if None in target_cds_positions:
            continue
        # sanity check - should probably remove. But should probably write tests too...
        assert all([target_cds_positions[2] - target_cds_positions[1] == 1, target_cds_positions[1] -
                    target_cds_positions[0] == 1, target_cds_positions[2] - target_cds_positions[0] == 2])
        target_codon = target_cds[target_cds_positions[0]:target_cds_positions[0] + 3]
        query_codon = query_cds[i:i + 3]
        assert len(target_codon) == len(query_codon) == 3, a.name
        yield target_cds_positions[0], target_codon, query_codon


def get_adjusted_starts_ends(t, aln):
    return [[aln.target_coordinate_to_query(intron.start - 1), aln.target_coordinate_to_query(intron.stop)]
            for intron in t.intron_intervals]


def is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance=5):
    q_gap_start = aln.target_coordinate_to_query(intron.start - 1)
    q_gap_end = aln.target_coordinate_to_query(intron.stop)
    return query_contains_intron(q_gap_start - fuzz_distance, q_gap_end + fuzz_distance, ref_starts)


def query_contains_intron(q_gap_start, q_gap_end, ref_starts):
    r = [q_gap_start <= ref_gap <= q_gap_end for ref_gap in ref_starts]
    return True if any(r) else False


def fix_ref_q_starts(ref_aln):
    if ref_aln.strand == "-":
        ref_starts = [ref_aln.q_size - (ref_aln.q_starts[i] + ref_aln.block_sizes[i]) for i in 
                      xrange(len(ref_aln.q_starts))]
    else:
        ref_starts = ref_aln.q_starts
    return ref_starts