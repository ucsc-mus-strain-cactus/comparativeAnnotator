"""
Extends align_augustus.py for consensus finding for comparative Augustus.
Takes every CDS for each transcript and maps against all transcripts for all genes assigned this transcript
(in in the name2 field)
Can be run in two modes - either aligning CGP transcripts (which are by definition CDS only) or extracting and aligning
the CDS of TM/TMR transcripts.
"""

import os
import argparse
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from lib.psl_lib import PslRow, remove_augustus_alignment_number, remove_alignment_number
from lib.sql_lib import ExclusiveSqlConnection, get_gene_transcript_map, attach_databases
from lib.seq_lib import GenePredTranscript
from lib.general_lib import tokenize_stream, grouper
from sonLib.bioio import fastaWrite, popenCatch, system, TempFileTree, catFiles
from pyfasta import Fasta
from lib.general_lib import format_ratio

__author__ = "Ian Fiddes"


def prepare_tmp_files(tmp_dir, gp, target_genome_fasta):
    """
    Builds temporary files for BLAT
    """
    tmp_tgt = os.path.join(tmp_dir, "tmp_cgp")
    tmp_ref = os.path.join(tmp_dir, "tmp_ref")
    tmp_psl = os.path.join(tmp_dir, "tmp_psl")
    cds = gp.get_cds(target_genome_fasta)
    with open(tmp_tgt, "w") as outf:
        outf.write(">{}\n{}\n".format(gp.name, cds))
    return tmp_tgt, tmp_ref, tmp_psl


def evaluate_blat_results(r):
    """
    Evalutes chained BLAT output for the one best alignment. Reports this alignments coverage and identity.
    """
    if len(r) == 0:
        return 0, 0
    else:
        p_list = [PslRow(x) for x in tokenize_stream(r)]
        # we take the smallest coverage value to account for Augustus adding bases
        p_list = [[min(x.coverage, x.target_coverage), x.identity] for x in p_list]
        best_cov, best_ident = sorted(p_list, key=lambda x: x[0])[-1]
        return best_cov, best_ident


def align_gp(target, genome, ref_genome, ref_tx_fasta, target_genome_fasta, gp, mode, out_db, comp_ann_path,
             chunk_size):
    """
    Initial wrapper job. Constructs a file tree and starts alignment job batches in groups of chunk_size.
    Follow on: concatenates file tree.
    """
    file_tree = TempFileTree(target.getGlobalTempDir())
    for recs in grouper(open(gp), chunk_size):
        target.addChildTargetFn(align_wrapper, args=[recs, file_tree, ref_tx_fasta, target_genome_fasta, comp_ann_path,
                                                     ref_genome, mode])
    target.setFollowOnTargetFn(cat, args=[genome, file_tree, out_db, mode])


def align_wrapper(target, recs, file_tree, ref_tx_fasta, target_genome_fasta, comp_ann_path, ref_genome, mode):
    """
    Alignment wrapper for grouped CGP records or grouped consensus records.
    For CGP mode, pulls down a gene -> transcript map and uses this to determine alignment targets, if they exist.
    """
    tmp_dir = target.getGlobalTempDir()
    results = []
    if mode == "cgp":
        con, cur = attach_databases(comp_ann_path, mode="reference")
        gene_transcript_map = get_gene_transcript_map(cur, ref_genome, biotype="protein_coding")
    for rec in recs:
        gp = GenePredTranscript(rec.rstrip().split("\t"))
        gene_names = gp.name2.split(",")
        if mode == "cgp":
            tx_dict = {n: gene_transcript_map[n] for n in gene_names if n in gene_transcript_map}
            if len(tx_dict) > 0:
                results.extend(align_cgp(tmp_dir, gp, target_genome_fasta, tx_dict, ref_tx_fasta))
        else:
            results.append(align_consensus(tmp_dir, gp, target_genome_fasta, ref_tx_fasta))
    with open(file_tree.getTempFile(), "w") as outf:
        for x in results:
            outf.write("".join([",".join(x), "\n"]))


def align_cgp(tmp_dir, gp, target_genome_fasta, tx_dict, ref_tx_fasta):
    """
    Main CGP alignment function. For each CGP transcript, uses tx_dict to BLAT against all transcripts. These alignments
    are then chained and the highest coverage alignment used. This circumvents problems with multiple self alignments
    in the case of repeats.
    """
    results = []
    ref_tx_fasta = Fasta(ref_tx_fasta)
    target_genome_fasta = Fasta(target_genome_fasta)
    tmp_tgt, tmp_ref, tmp_psl = prepare_tmp_files(tmp_dir, gp, target_genome_fasta)
    for gene_name, tx_names in tx_dict.iteritems():
        for tx_name in tx_names:
            tx_seq = str(ref_tx_fasta[tx_name])
            fastaWrite(tmp_ref, tx_name, tx_seq)
            system("blat {} {} -out=psl -noHead {}".format(tmp_tgt, tmp_ref, tmp_psl))
            r = popenCatch("simpleChain -outPsl {} /dev/stdout".format(tmp_psl))
            r = r.split("\n")[:-1]
            best_cov, best_ident = evaluate_blat_results(r)
            results.append(map(str, [gp.name, gene_name, tx_name, best_cov, best_ident]))
    return results


def align_consensus(tmp_dir, gp, target_genome_fasta, ref_tx_fasta):
    """
    Main consensus alignment function.
    """
    ref_tx_fasta = Fasta(ref_tx_fasta)
    target_genome_fasta = Fasta(target_genome_fasta)
    tmp_tgt, tmp_ref, tmp_psl = prepare_tmp_files(tmp_dir, gp, target_genome_fasta)
    tx_seq = str(ref_tx_fasta[gp.name])
    fastaWrite(tmp_ref, gp.name, tx_seq)
    system("blat {} {} -out=psl -noHead {}".format(tmp_tgt, tmp_ref, tmp_psl))
    r = popenCatch("simpleChain -outPsl {} /dev/stdout".format(tmp_psl))
    r = r.split("\n")[:-1]
    best_cov, best_ident = evaluate_blat_results(r)
    return map(str, [gp.id, gp.name, best_cov, best_ident])


def cat(target, genome, file_tree, out_db, mode):
    """
    Concatenates the final resulting file tree before database construction.
    """
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[genome, tmp_file, out_db, mode])


def load_db(target, genome, tmp_file, out_db, mode):
    """
    Loads the data into a sqlite database.
    """
    if mode == "cgp":
        names = ["CgpId", "GeneId", "EnsId", "AlignmentCoverage", "AlignmentIdentity"]
        df = pd.read_csv(tmp_file, index_col=[0, 2], names=names)
    else:
        names = ["AlignmentId", "EnsId", "AlignmentCoverage", "AlignmentIdentity"]
        df = pd.read_csv(tmp_file, index_col=0, names=names)
    df = df.convert_objects(convert_numeric=True)  # have to convert to float because pandas lacks a good dtype function
    df = df.sort_index()
    table = "_".join([genome, mode])
    with ExclusiveSqlConnection(out_db) as con:
        df.to_sql(table, con, if_exists="replace", index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refGenome", required=True)
    parser.add_argument("--refTranscriptFasta", required=True)
    parser.add_argument("--targetGenomeFasta", required=True)
    parser.add_argument("--outDb", default="cgp_cds_metrics.db")
    parser.add_argument("--compAnnPath", required=True)
    gp_group = parser.add_mutually_exclusive_group(required=True)
    gp_group.add_argument("--cgpGp")
    gp_group.add_argument("--consensusGp")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    out_db = os.path.join(args.compAnnPath, args.outDb)
    if args.cgpGp is not None:
        gp = args.cgpGp
        mode = "cgp"
        chunk_size = 15  # smaller chunk size because we will do more alignments per transcript
    else:
        gp = args.consensusGp
        mode = "consensus"
        chunk_size = 40
    s = Stack(Target.makeTargetFn(align_gp, args=[args.genome, args.refGenome, args.refTranscriptFasta, 
                                                  args.targetGenomeFasta, gp, mode, out_db, args.compAnnPath,
                                                  chunk_size]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from scripts.align_cgp_cds import *
    main()
