import os
import argparse
import itertools
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from lib.seq_lib import get_sequence_dict
from lib.psl_lib import PslRow, remove_augustus_alignment_number, remove_alignment_number
from lib.general_lib import tokenize_stream, format_ratio
from lib.sql_lib import ExclusiveSqlConnection
from sonLib.bioio import fastaWrite, popenCatch, TempFileTree, catFiles


def coverage(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.mismatches for x in p_list)
    rep = sum(x.repmatches for x in p_list)
    # ident/cov can end up slightly above 1 due to adding floats
    cov = 100 * format_ratio(m + mi + rep, p_list[0].q_size)
    return min(cov, 100.0)


def identity(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.mismatches for x in p_list)
    rep = sum(x.repmatches for x in p_list)
    ins = sum(x.q_num_insert for x in p_list)
    # ident/cov can end up slightly above 1 due to adding floats
    ident = 100 * format_ratio(m + rep, m + rep + mi + ins)
    return min(ident, 100.0)


def seq_dict_paired_chunker(ref_dict, tx_dict, size=200):
    target_it = iter(tx_dict)
    for i in xrange(0, len(tx_dict), size):
        yield [[tx_dict[aln_id], ref_dict[psl_lib.strip_alignment_numbers(aln_id)]] for aln_id in 
               itertools.islice(target_it, size)]


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
            results.append([aug_aln_id, aln_id, "0", "0"])
        else:
            p_list = [PslRow(x) for x in tokenize_stream(r)]
            results.append(map(str, [aug_aln_id, aln_id, identity(p_list), coverage(p_list)]))
    with open(file_tree.getTempFile(), "w") as outf:
        for x in results:
            outf.write("".join([",".join(x), "\n"]))


def align(target, file_tree, chunk, ref_fasta, target_fasta):
    ref_seq_dict = get_sequence_dict(ref_fasta)
    target_seq_dict = get_sequence_dict(target_fasta)
    results = []
    for tx, ref_tx in chunk:
        tmp_results = [tx.name, remove_augustus_alignment_number(tx.name)]
        for mode in ["mrna", "cds"]:
            tmp_aug = os.path.join(target.getLocalTempDir(), "tmp_aug_{}".format(mode))
            tmp_gencode = os.path.join(target.getLocalTempDir(), "tmp_gencode_{}".format(mode))
            tx_fn = eval("tx.get_{}".format(mode))
            tx_seq = tx_fn(target_seq_dict)
            ref_fn = eval("ref_tx.get_{}".format(mode))
            ref_seq = ref_fn(ref_seq_dict)
            fastaWrite(tmp_aug, tx.name, tx_seq)
            fastaWrite(tmp_gencode, tx.name, ref_seq)
            r = popenCatch("blat {} {} -out=psl -noHead /dev/stdout".format(tmp_gencode, tmp_aug))
            r = r.split("\n")[:-3]
            if len(r) == 0:
                tmp_results.extend(["0", "0"])
            else:
                p_list = [PslRow(x) for x in tokenize_stream(r)]
                tmp_results.extend([identity(p_list), coverage(p_list)])
        results.append(tmp_results)
    with open(file_tree.getTempFile(), "w") as outf:
        for x in results:
            outf.write("".join([",".join(x), "\n"]))


def align_augustus(target, genome, ref_fasta, target_fasta, ref_gp, aug_gp, out_dir):
    file_tree = TempFileTree(target.getGlobalTempDir()))
    ref_dict = seq_lib.get_transcript_dict(ref_gp)
    tx_dict = seq_lib.get_transcript_dict(aug_gp)
    for chunk in seq_dict_paired_chunker(ref_dict, tx_dict, size=200):
        target.addChildTargetFn(align, args=[file_tree, chunk, ref_fasta, target_fasta])
    target.setFollowOnTargetFn(cat, args=(genome, file_tree, out_dir))


def cat(target, genome, file_tree, out_dir):
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[genome, tmp_file, out_dir])


def load_db(target, genome, tmp_file, out_dir):
    df = pd.read_csv(tmp_file, index_col=0, names=["AugustusAlignmentId", "AlignmentId", "AlignmentIdentity", 
                                                   "AlignmentCoverage", "CdsAlignmentIdentity", "CdsAlignmentCoverage"])
    df = df.convert_objects(convert_numeric=True)  # have to convert to float because pandas lacks a good dtype function
    df = df.sort_index()
    database_path = os.path.join(out_dir, "augustus_attributes.db")
    with ExclusiveSqlConnection(database_path) as con:
        df.to_sql(genome, con, if_exists="replace", index_label="AugustusAlignmentId")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refFasta", required=True)
    parser.add_argument("--targetFasta", required=True)
    parser.add_argument("--refGp", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--outDir", required=True)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    i = Stack(Target.makeTargetFn(align_augustus, args=[args.genome, args.refFasta, args.targetFasta, args.refGp,
                                                        args.augGp, args.outDir])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.align_augustus import *
    main()
