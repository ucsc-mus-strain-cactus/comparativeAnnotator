import sys
import os
import argparse
import itertools
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.comp_ann_lib as comp_ann_lib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, help="genome to get intron information from")
    parser.add_argument("--psl", required=True, help="transMap PSL")
    parser.add_argument("--refPsl", required=True, help="reference fake PSL")
    parser.add_argument("--gp", required=True)
    parser.add_argument("--outPath", required=True, help="File name for output")
    parser.add_argument("--fuzz_distance", default=5, help="how much fuzz_distance do we allow?")
    return parser.parse_args()


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


def main():
    args = parse_args()
    aln_dict = psl_lib.get_alignment_dict(args.psl)
    ref_aln_dict = psl_lib.get_alignment_dict(args.refPsl)
    tx_dict = seq_lib.get_transcript_dict(args.gp)
    with open(args.outPath, "w") as outf:
        for aln_id, aln in sorted(aln_dict.iteritems(), key=lambda x: x[0]):
            ref_aln = ref_aln_dict[psl_lib.remove_alignment_number(aln_id)]
            t = tx_dict[aln_id]
            vec = build_intron_vector(aln, ref_aln, t, args.fuzz_distance)
            outf.write("{}\t{}\n".format(aln_id, ",".join(vec)))


if __name__ == "__main__":
    main()
