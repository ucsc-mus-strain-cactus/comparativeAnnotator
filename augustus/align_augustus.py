"""
Aligns AugustusTMR transcripts to their respective reference transMap transcripts, producing coverage and identity
metrics in a sqlite database. This is used for building consensus gene sets.
"""

import os
import pandas as pd
from collections import defaultdict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from pycbio.bio.psl import PslRow
from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number, remove_augustus_alignment_number, strip_alignment_numbers
from pycbio.sys.dataOps import grouper
from pycbio.sys.fileOps import iterRows, ensureDir
from pycbio.bio.transcripts import get_transcript_dict
from sonLib.bioio import fastaWrite, popenCatch, system, TempFileTree, catFiles
from pyfasta import Fasta
from pycbio.sys.sqliteOps import ExclusiveSqlConnection


def find_best_overall_aln(results, cov_weight=0.25, ident_weight=0.75):
    weighted_results = []
    for tgt_id, query_id, gencode_id, best_cov, best_ident, num_introns in results:
        weighted_results.append([tgt_id, cov_weight * best_cov + ident_weight * best_ident])
    best = sorted(weighted_results, key=lambda (aln_id, s): -s)[0][0]
    overall_results = []
    for tgt_id, query_id, gencode_id, best_cov, best_ident, num_introns in results:
        is_best = 1 if tgt_id == best else 0
        overall_results.append([tgt_id, query_id, gencode_id, best_cov, best_ident, num_introns, is_best])
    return overall_results


def align(target, target_fasta, chunk, ref_fasta, file_tree):
    g_f = Fasta(target_fasta)
    r_f = Fasta(ref_fasta)
    results = []
    tmp_aug = os.path.join(target.getGlobalTempDir(), "tmp_aug")
    tmp_gencode = os.path.join(target.getGlobalTempDir(), "tmp_gencode")
    tmp_psl = os.path.join(target.getGlobalTempDir(), 'tmp_psl')
    for tx_id, vals in chunk:
        these_results = []
        for tgt_id, tx in vals:
            query_id = remove_augustus_alignment_number(tgt_id)
            gencode_id = remove_alignment_number(query_id)
            gencode_seq = str(r_f[gencode_id])
            aug_seq = str(g_f[tgt_id])
            fastaWrite(tmp_aug, tgt_id, aug_seq)
            fastaWrite(tmp_gencode, gencode_id, gencode_seq)
            system("blat {} {} -out=psl -noHead {}".format(tmp_aug, tmp_gencode, tmp_psl))
            r = popenCatch("simpleChain -outPsl {} /dev/stdout".format(tmp_psl))
            r = r.split("\n")[:-1]
            if len(r) == 0:
                these_results.append([tgt_id, query_id, gencode_id, 0, 0, len(tx.intron_intervals)])
            else:
                p_list = [PslRow(x) for x in iterRows(r)]
                # we take the smallest coverage value to account for Augustus adding bases
                p_list = [[min(x.coverage, x.target_coverage), x.identity] for x in p_list]
                best_cov, best_ident = sorted(p_list, key=lambda x: x[0])[-1]
                these_results.append([tgt_id, query_id, gencode_id, best_cov, best_ident, len(tx.intron_intervals)])
        results.extend(find_best_overall_aln(these_results))
    with open(file_tree.getTempFile(), "w") as outf:
        for x in results:
            outf.write("".join([",".join(map(str, x)), "\n"]))


def get_hierarchical_transcript_dict(gp):
    tx_dict = get_transcript_dict(gp)
    h_dict = defaultdict(list)
    for tx_id, tx in tx_dict.iteritems():
        h_dict[strip_alignment_numbers(tx_id)].append([tx_id, tx])
    return h_dict


def align_augustus_wrapper(target, args):
    file_tree = TempFileTree(target.getGlobalTempDir())
    h_dict = get_hierarchical_transcript_dict(args.augustus_gp)
    for chunk in grouper(h_dict.iteritems(), 50):
        target.addChildTargetFn(align, args=[args.fasta, chunk, args.ref_fasta, file_tree])
    target.setFollowOnTargetFn(cat, args=(args.genome, file_tree, args.db))


def cat(target, genome, file_tree, out_db):
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[genome, tmp_file, out_db])


def load_db(target, genome, tmp_file, out_db):
    df = pd.read_csv(tmp_file, index_col=0, names=["AugustusAlignmentId", "AlignmentId", "TranscriptId",
                                                   "AugustusAlignmentCoverage", "AugustusAlignmentIdentity",
                                                   "AugustusNumberIntrons", "AugustusBestOverallAln"])
    df = df.convert_objects(convert_numeric=True)  # have to convert to float because pandas lacks a good dtype function
    df = df.sort_index()
    ensureDir(os.path.dirname(out_db))
    with ExclusiveSqlConnection(out_db) as con:
        df.to_sql(genome + '_AugustusAttributes', con, if_exists="replace", index_label='AugustusAlignmentId')


def align_augustus(args):
    i = Stack(Target.makeTargetFn(align_augustus_wrapper, args=[args])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
