"""
Simple script that reports 1 for introns that existed in the original alignment, and 0 for those that did not.
"""
import comparativeAnnotator.comp_lib.annotation_utils as comp_ann_lib
from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.database_queries import get_transcript_biotype_map


def build_intron_vector(aln, ref_aln, t, fuzz_distance):
    result = []
    ref_starts = comp_ann_lib.fix_ref_q_starts(ref_aln)
    for intron in t.intron_intervals:
        if comp_ann_lib.short_intron(intron):
            result.append("0")
        elif comp_ann_lib.is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance) is False:
            result.append("0")
        else:
            result.append("1")
    return result


def find_intron_vector(args, fuzz_distance=10):
    aln_dict = get_alignment_dict(args.psl)
    ref_aln_dict = get_alignment_dict(args.ref_psl)
    tx_dict = get_transcript_dict(args.gp)
    tx_biotype_map = get_transcript_biotype_map(args.query_genome, args.db)
    for aln_id, aln in sorted(aln_dict.iteritems(), key=lambda x: x[0]):
        if tx_biotype_map[remove_alignment_number(aln_id)] != 'protein_coding':
            continue
        ref_aln = ref_aln_dict[remove_alignment_number(aln_id)]
        t = tx_dict[aln_id]
        vec = build_intron_vector(aln, ref_aln, t, fuzz_distance)
        yield '\t'.join(map(str, t.get_gene_pred()) + [','.join(vec)]) + '\n'
