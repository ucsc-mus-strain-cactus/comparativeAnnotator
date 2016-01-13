"""
Aligns AugustusTMR transcripts to their respective reference transMap transcripts, producing coverage and identity
metrics in a sqlite database. This is used for building consensus gene sets.
"""

import os
import argparse
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from lib.psl_lib import PslRow, remove_augustus_alignment_number, remove_alignment_number
from lib.general_lib import tokenize_stream, grouper
from sonLib.bioio import fastaWrite, popenCatch, system, TempFileTree, catFiles
from pyfasta import Fasta
from lib.general_lib import format_ratio
from lib.sql_lib import ExclusiveSqlConnection


def align(target, target_fasta, chunk, ref_fasta, file_tree):
    g_f = Fasta(target_fasta)
    r_f = Fasta(ref_fasta)
    results = []
    for tgt_id in chunk:
        query_id = remove_augustus_alignment_number(tgt_id)
        gencode_id = remove_alignment_number(query_id)
        gencode_seq = str(r_f[gencode_id])
        aug_seq = str(g_f[tgt_id])
        tmp_aug = os.path.join(target.getLocalTempDir(), "tmp_aug")
        tmp_gencode = os.path.join(target.getLocalTempDir(), "tmp_gencode")
        fastaWrite(tmp_aug, tgt_id, aug_seq)
        fastaWrite(tmp_gencode, gencode_id, gencode_seq)
        tmp_psl = os.path.join(target.getLocalTempDir(), tgt_id)
        system("blat {} {} -out=psl -noHead {}".format(tmp_aug, tmp_gencode, tmp_psl))
        r = popenCatch("simpleChain -outPsl {} /dev/stdout".format(tmp_psl))
        r = r.split("\n")[:-1]
        if len(r) == 0:
            results.append([tgt_id, query_id, "0", "0"])
        else:
            p_list = [PslRow(x) for x in tokenize_stream(r)]
            # we take the smallest coverage value to account for Augustus adding bases
            p_list = [[min(x.coverage, x.target_coverage), x.identity] for x in p_list]
            best_cov, best_ident = sorted(p_list, key=lambda x: x[0])[-1]
            results.append(map(str, [tgt_id, query_id, best_cov, best_ident]))
    with open(file_tree.getTempFile(), "w") as outf:
        for x in results:
            outf.write("".join([",".join(x), "\n"]))


def align_augustus(target, genome, ref_fasta, target_fasta, target_fasta_index, out_db):
    file_tree = TempFileTree(target.getGlobalTempDir())
    tgt_ids = [x.split()[0] for x in open(target_fasta_index)]
    for chunk in grouper(tgt_ids, 200):
        target.addChildTargetFn(align, args=[target_fasta, chunk, ref_fasta, file_tree])
    target.setFollowOnTargetFn(cat, args=(genome, file_tree, out_db))


def cat(target, genome, file_tree, out_db):
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[genome, tmp_file, out_db])


def load_db(target, genome, tmp_file, out_db):
    df = pd.read_csv(tmp_file, index_col=0, names=["AugustusAlignmentId", "AlignmentId", "AlignmentIdentity", 
                                                   "AlignmentCoverage"])
    df = df.convert_objects(convert_numeric=True)  # have to convert to float because pandas lacks a good dtype function
    df = df.sort_index()
    with ExclusiveSqlConnection(out_db) as con:
        df.to_sql(genome, con, if_exists="replace", index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refTranscriptFasta", required=True)
    parser.add_argument("--targetTranscriptFasta", required=True)
    parser.add_argument("--targetTranscriptFastaIndex", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--outDb", default="augustus_attributes.db")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    out_db = os.path.join(args.outDir, args.outDb)
    i = Stack(Target.makeTargetFn(align_augustus, args=[args.genome, args.refTranscriptFasta, 
                                                        args.targetTranscriptFasta, args.targetTranscriptFastaIndex,
                                                        out_db])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.align_augustus import *
    main()