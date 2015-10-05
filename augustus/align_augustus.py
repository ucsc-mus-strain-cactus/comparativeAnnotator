import os
import argparse
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from lib.psl_lib import PslRow, remove_augustus_alignment_number, remove_alignment_number
from sonLib.bioio import fastaWrite, popenCatch, TempFileTree, catFiles
from pyfaidx import Fasta
from lib.general_lib import format_ratio
from lib.sql_lib import ExclusiveSqlConnection


def coverage(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.mismatches for x in p_list)
    rep = sum(x.repmatches for x in p_list)
    # ident/cov can end up slightly above 1 due to adding floats
    cov = format_ratio(m + mi + rep, p_list[0].q_size)
    return min(cov, 1.0)


def identity(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.mismatches for x in p_list)
    rep = sum(x.repmatches for x in p_list)
    ins = sum(x.q_num_insert for x in p_list)
    # ident/cov can end up slightly above 1 due to adding floats
    ident = format_ratio(m + rep, m + rep + mi + ins)
    return min(ident, 1.0)


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def align(target, target_fasta, chunk, ref_fasta, file_tree):
    g_f = Fasta(target_fasta)
    r_f = Fasta(ref_fasta)
    results = []
    for aug_aln_id in chunk:
        aln_id = remove_augustus_alignment_number(aug_aln_id)
        gencode_id = remove_alignment_number(aln_id)
        gencode_seq = str(r_f[gencode_id])
        aug_seq = str(g_f[aug_aln_id])
        tmp_aug = os.path.join(target.getLocalTempDir(), "tmp_aug")
        tmp_gencode = os.path.join(target.getLocalTempDir(), "tmp_gencode")
        fastaWrite(tmp_aug, aug_aln_id, aug_seq)
        fastaWrite(tmp_gencode, gencode_id, gencode_seq)
        r = popenCatch("blat {} {} -out=psl -noHead /dev/stdout".format(tmp_gencode, tmp_aug))
        r = r.split("\n")[:-3]
        if len(r) == 0:
            results.append([aug_aln_id, "0", "0"])
        else:
            p_list = [PslRow(x) for x in r]
            results.append(map(str, [aug_aln_id, identity(p_list), coverage(p_list)]))
    with open(file_tree.getTempFile(), "w") as outf:
        outf.write("AlignmentId,AlignmentIdentity,AlignmentCoverage\n")
        for x in results:
            outf.write("".join([",".join(x), "\n"]))


def align_augustus(target, genome, ref_fasta, target_fasta, target_fasta_index, out_dir):
    file_tree = TempFileTree(target.getGlobalTempDir())
    aug_aln_ids = [x.split()[0] for x in open(target_fasta_index)]
    for chunk in chunker(aug_aln_ids, 200):
        target.addChildTargetFn(align, args=[target_fasta, chunk, ref_fasta, file_tree])
    target.setFollowOnTargetFn(cat, args=(genome, file_tree, out_dir))


def cat(target, genome, file_tree, out_dir):
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[genome, tmp_file, out_dir])


def load_db(target, genome, tmp_file, out_dir):
    df = pd.DataFrame.from_csv(tmp_file)
    database_path = os.path.join(out_dir, "augustus_attributes.db")
    df = df.sort_index()
    with ExclusiveSqlConnection(database_path) as con:
        df.to_sql(genome, con, if_exists="replace", index_label="AlignmentId")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refTranscriptFasta", required=True)
    parser.add_argument("--targetTranscriptFasta", required=True)
    parser.add_argument("--targetTranscriptFastaIndex", required=True)
    parser.add_argument("--outDir", required=True)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    i = Stack(Target.makeTargetFn(align_augustus, args=[args.genome, args.refTranscriptFasta, 
                                                        args.targetTranscriptFasta, args.targetTranscriptFastaIndex,
                                                        args.outDir])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.align_augustus import *
    main()
