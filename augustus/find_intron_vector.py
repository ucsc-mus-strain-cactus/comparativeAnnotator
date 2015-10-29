import sys
import os
import argparse
import itertools
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, help="genome to get intron information from")
    parser.add_argument("--gp", required=True, help="genePred for this genome's transMap results")
    parser.add_argument("--psl", required=True, help="transMap PSL")
    parser.add_argument("--refGp", required=True, help="reference genePred")
    parser.add_argument("--outPath", required=True, help="File name for output")
    return parser.parse_args()

"""
def build_intron_vector(a, t, aln):
    result = []
    offset = aln.target_coordinate_to_query(t.transcript_coordinate_to_chromosome(0))
    target_splices = {x.stop - offset for x in t.exons[:-1]}
    query_splices = {x.stop for x in a.exons[:-1]}
    not_original_splices = target_splices - query_splices
    for exon in t.exons[:-1]:
        target_splice = exon.stop + offset
        if target_splice in not_original_splices:
            result.append(0)
        else:
            result.append(1)
    return result
"""
def build_intron_vector(a, t, aln):
    result = []
    for intron in t.intron_intervals:
        if len(intron) <= 30:
            result.append(0)
        else:
            result.append(1)


def main():
    args = parse_args()
    ref_dict = seq_lib.get_transcript_dict(args.refGp)
    target_dict = seq_lib.get_transcript_dict(args.gp)
    aln_dict = psl_lib.get_alignment_dict(args.psl)
    with open(args.outPath, "w") as outf:
        for aln_id, t in sorted(target_dict.iteritems(), key=lambda x: x[0]):
            a = ref_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = aln_dict[aln_id]
            vec = build_intron_vector(a, t, aln)
            outf.write("{}\t{}\n".format(aln_id, ",".join(map(str, vec))))


if __name__ == "__main__":
    main()
