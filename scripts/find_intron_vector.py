import sys
import os
import argparse
from scripts.plot_functions import attach_databases
import lib.sequence_lib as seq_lib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, help="genome to get intron information from")
    parser.add_argument("--gp", required=True, help="genePred for this genome's transMap results")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--outPath", required=True, help="File name for output.")
    return parser.parse_args()


def load_gp(path):
    tm_recs = seq_lib.getGenePredTranscripts(path)
    tm_dict = seq_lib.transcriptListToDict(tm_recs, noDuplicates=True)
    return tm_dict


def load_database_results(cur, genome):
    cmd = "SELECT AlignmentId, HasOriginalIntrons from details.'{}'".format(genome)
    results = {x: y for x, y in cur.execute(cmd).fetchall()}
    return results


def parse_db_rec(db_rec):
    if db_rec is not None:
        db_rec = db_rec.split("\n")
        result = {tuple(map(int, x.split()[1:3])) for x in db_rec}
        return result
    return set()


def build_intron_vector(gene_rec, db_rec):
    db_rec_set = parse_db_rec(db_rec)
    result = []
    for intron in gene_rec.intronIntervals:
        if len(intron) == 0:
            # ignore the 0bp introns that Mark has in the new chaining
            continue
        if (intron.start, intron.stop) in db_rec_set:
            result.append(0)
        else:
            result.append(1)
    return result


def main():
    args = parse_args()
    con, cur = attach_databases(args.comparativeAnnotationDir)
    tm_dict = load_gp(args.gp)
    db_dict = load_database_results(cur, args.genome)
    with open(args.outPath, "w") as outf:
        for aln_id, gene_rec in sorted(tm_dict.iteritems(), key=lambda x: x[0]):
            db_rec = db_dict[aln_id]
            vec = build_intron_vector(gene_rec, db_rec)
            outf.write("{}\t{}\n".format(aln_id, ",".join(map(str, vec))))


if __name__ == "__main__":
    main()
