import sys
import os
import argparse
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


def load_gp(path):
    tm_recs = seq_lib.get_gene_pred_transcripts(path)
    tm_dict = seq_lib.transcript_list_to_dict(tm_recs)
    return tm_dict


def load_psl(path):
    psl_recs = psl_lib.read_psl(path)
    psl_dict = psl_lib.get_psl_dict(psl_recs)
    return psl_dict


def build_intron_vector(a, t, aln):
    original_introns = {(x.start, x.stop) for x in a.intron_intervals}
    result = []
    for intron in t.intron_intervals:
        a_start = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.start - 1)) + 1
        a_stop = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.stop))
        if (a.start, a.stop) not in original_introns:
            result.append(0)
        else:
            result.append(1)
    return result


def main():
    args = parse_args()
    ref_dict = load_gp(args.refGp)
    target_dict = load_gp(args.gp)
    aln_dict = load_psl(args.psl)
    with open(args.outPath, "w") as outf:
        for aln_id, t in sorted(target_dict.iteritems(), key=lambda x: x[0]):
            a = ref_dict[psl_lib.remove_alignment_number(aln_id)]
            aln = aln_dict[aln_id]
            vec = build_intron_vector(a, t, aln)
            outf.write("{}\t{}\n".format(aln_id, ",".join(map(str, vec))))


if __name__ == "__main__":
    main()
