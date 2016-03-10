"""
Extends align_augustus.py for consensus finding for comparative Augustus.
Takes every CDS for each transcript and maps against all transcripts for all genes assigned this transcript
(in in the name2 field)
Can be run in two modes - either aligning CGP transcripts (which are by definition CDS only) or extracting and aligning
the CDS of TM/TMR transcripts.
"""

import os
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from pycbio.bio.psl import PslRow
from pycbio.bio.transcripts import GenePredTranscript
from pycbio.sys.sqliteOps import ExclusiveSqlConnection, execute_query
from pycbio.sys.fileOps import iterRows, ensureFileDir
from pycbio.sys.dataOps import grouper
from comparativeAnnotator.database_queries import get_gene_transcript_map
from sonLib.bioio import fastaWrite, popenCatch, system, TempFileTree, catFiles
from pyfasta import Fasta

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
        p_list = [PslRow(x) for x in iterRows(r)]
        # we take the smallest coverage value to account for Augustus adding bases
        p_list = [[min(x.coverage, x.target_coverage), x.identity] for x in p_list]
        best_cov, best_ident = sorted(p_list, key=lambda x: x[0])[-1]
        return best_cov, best_ident


def align_wrapper(target, recs, file_tree, ref_tx_fasta, target_genome_fasta, comp_db_path, ref_genome, mode):
    """
    Alignment wrapper for grouped CGP records or grouped consensus records.
    For CGP mode, pulls down a gene -> transcript map and uses this to determine alignment targets, if they exist.
    """
    tmp_dir = target.getGlobalTempDir()
    results = []
    if mode == "cgp":
        gene_transcript_map = get_gene_transcript_map(ref_genome, comp_db_path, biotype="protein_coding")
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


def cat(target, file_tree, out_db, mode, table):
    """
    Concatenates the final resulting file tree before database construction.
    """
    tmp_file = os.path.join(target.getGlobalTempDir(), "tmp.txt")
    catFiles(file_tree.listFiles(), tmp_file)
    target.setFollowOnTargetFn(load_db, args=[tmp_file, out_db, mode, table])


def load_db(target, tmp_file, out_db, mode, table):
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
    ensureFileDir(out_db)
    with ExclusiveSqlConnection(out_db) as con:
        df.to_sql(table, con, if_exists="replace", index=True)
        cur = con.cursor()
        cmd = 'CREATE TABLE IF NOT EXISTS completionFlags (aln_table TEXT)'
        execute_query(cur, cmd)
        cmd = 'INSERT INTO completionFlags (aln_table) VALUES ("{}")'.format(table)
        execute_query(cur, cmd)


def align_gp(target, ref_genome, ref_tx_fasta, target_genome_fasta, gp, mode, out_db, comp_db_path,
             table, chunk_size):
    """
    Initial wrapper job. Constructs a file tree and starts alignment job batches in groups of chunk_size.
    Follow on: concatenates file tree.
    """
    file_tree = TempFileTree(target.getGlobalTempDir())
    for recs in grouper(open(gp), chunk_size):
        target.addChildTargetFn(align_wrapper, args=[recs, file_tree, ref_tx_fasta, target_genome_fasta, comp_db_path,
                                                     ref_genome, mode])
    target.setFollowOnTargetFn(cat, args=[file_tree, out_db, mode, table])


def align_cgp_cds(args):
    assert args.mode in ['cgp', 'consensus']
    chunk_size = 20 if args.mode == 'cgp' else 75
    s = Stack(Target.makeTargetFn(align_gp, args=[args.refGenome, args.refTranscriptFasta,
                                                  args.targetGenomeFasta, args.gp, args.mode, args.cgpDb,
                                                  args.compDb, args.table, chunk_size]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
