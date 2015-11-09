"""
Script to build the databases and tracks from comparativeAnnotator results
"""
import os
import subprocess
import cPickle as pickle
import pandas as pd

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import mkdir_p
import etc.config

__author__ = "Ian Fiddes"

augustus_classifier_tracks = [etc.config.allAugustusClassifiers]
ref_tracks = [etc.config.refClassifiers]
tm_classifier_tracks = [etc.config.allClassifiers, etc.config.potentiallyInterestingBiology, etc.config.assemblyErrors,
                        etc.config.alignmentErrors]


def database_wrapper(target, args, tmp_dir):
    """
    Calls database for each database in this analysis.
    """
    if args.mode == "augustus":
        for db in ["classify", "details"]:
            db_path = os.path.join(args.outDir, "augustus_{}.db".format(db))
            database(args.genome, db, db_path, tmp_dir, args.mode)
    elif args.mode == "reference":
        for db in ["classify", "details"]:
            db_path = os.path.join(args.outDir, "{}.db".format(db))
            database(args.refGenome, db, db_path, tmp_dir, args.mode)
        attr_db_path = os.path.join(args.outDir, "attributes.db")
        ref_attr_table(args.refGenome, attr_db_path, args.gencodeAttributes, args.annotationGp)
    elif args.mode == "transMap":
        for db in ["classify", "details", "attributes"]:
            db_path = os.path.join(args.outDir, "{}.db".format(db))
            database(args.genome, db, db_path, tmp_dir, args.mode)
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")
    target.setFollowOnTargetFn(build_tracks_wrapper, args=[args])


def database(genome, db, db_path, tmp_dir, mode):
    data_dict = {}
    mkdir_p(os.path.dirname(db_path))
    data_path = os.path.join(tmp_dir, db)
    for col in os.listdir(data_path):
        p = os.path.join(data_path, col)
        with open(p) as p_h:
            data_dict[col] = pickle.load(p_h)
    if mode == "reference":
        index_label = "TranscriptId"
    elif mode == "transMap":
        index_label = "AlignmentId"
    else:
        index_label = "AugustusAlignmentId"
        # Hack to add transMap alignment ID column to Augustus databases.
        aug_ids = data_dict.itervalues().next().viewkeys()
        data_dict["AlignmentId"] = {x: psl_lib.remove_augustus_alignment_number(x) for x in aug_ids}
    sql_lib.write_dict(data_dict, db_path, genome, index_label)


def ref_attr_table(ref_genome, db_path, attr_file, ref_gp):
    """
    This function is used to add an extra table in reference mode holding all of the basic attributes.
    Basically directly dumping the tsv into sqlite3 with the addition of a refChrom column.
    """
    df = pd.read_table(attr_file, sep="\t", index_col=3, header=0)
    ref_dict = seq_lib.get_transcript_dict(ref_gp)
    chromosome_dict = {"refChrom": {x: y.chromosome for x, y in ref_dict.iteritems()}}
    chromosome_df = pd.DataFrame.from_dict(chromosome_dict)
    df2 = pd.merge(df, chromosome_df, left_index=True, right_index=True)
    with sql_lib.ExclusiveSqlConnection(db_path) as con:
        df2.to_sql(ref_genome, con, if_exists="replace", index_label="TranscriptId")


def build_tracks_wrapper(target, args):
    if args.mode == "reference":
        classifier_tracks = ref_tracks
        genome = args.refGenome
    elif args.mode == "augustus":
        classifier_tracks = augustus_classifier_tracks
        genome = args.genome
    elif args.mode == "transMap":
        classifier_tracks = tm_classifier_tracks
        genome = args.genome
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")
    for query_fn in classifier_tracks:
        target.addChildTargetFn(build_classifier_tracks, args=[query_fn, genome, args])
    target.addChildTargetFn(build_good_track, args=[args])


def get_bed_paths(out_dir, query_name, genome):
    out_bed_dir = os.path.join(out_dir, "bedfiles", query_name, genome)
    out_bed_path = os.path.join(out_bed_dir, "{}.bed".format(genome))
    out_big_bed_dir = os.path.join(out_dir, "bigBedfiles", query_name, genome)
    out_big_bed_path = os.path.join(out_big_bed_dir, "{}.bb".format(genome))
    mkdir_p(out_bed_dir)
    mkdir_p(out_big_bed_dir)
    return out_bed_path, out_big_bed_path


def make_big_bed(out_bed_path, sizes, out_big_bed_path):
    subprocess.call(["bedSort", out_bed_path, out_bed_path])
    subprocess.call(["bedToBigBed", "-extraIndex=name", out_bed_path, sizes, out_big_bed_path])


def build_classifier_tracks(target, query_fn, genome, args):
    query = query_fn(genome)
    query_name = query_fn.__name__
    con, cur = sql_lib.attach_databases(args.outDir, mode=args.mode)
    bed_recs = cur.execute(query)
    out_bed_path, out_big_bed_path = get_bed_paths(args.outDir, query_name, genome)
    with open(out_bed_path, "w") as outf:
        for recs in bed_recs:
            for rec in recs:
                if rec is not None:
                    outf.write(rec)
    make_big_bed(out_bed_path, args.sizes, out_big_bed_path)


def get_all_tm_good(cur, ref_genome, genome):
    """
    transMap Good varies depending on if the transcript is coding or noncoding. We will build a set of IDs for both.
    """
    biotypes = sql_lib.get_all_biotypes(cur, genome, gene_level=False)
    all_ids = set()
    for biotype in biotypes:
        query = etc.config.transMapEval(ref_genome, genome, biotype=biotype, good=False)
        query_ids = sql_lib.get_query_ids(cur, query)
        all_ids |= query_ids
    return all_ids


def build_good_track(target, args):
    """
    Builds a specific track of Good transcripts for the current mode.
    """
    colors = {"coding": "59,101,69", "noncoding": "98,124,191", "not_good": "152,0,67"}
    con, cur = sql_lib.attach_databases(args.outDir, args.mode)
    biotype_map = sql_lib.get_transcript_biotype_map(cur, args.refGenome)
    if args.mode == "augustus":
        query = etc.config.augustusEval(args.genome, args.refGenome)
        good_ids = sql_lib.get_query_ids(cur, query)
        out_good_bed_path, out_good_big_bed_path = get_bed_paths(args.outDir, "augustus", args.genome)
        gp_dict = seq_lib.get_transcript_dict(args.augustusGp)
    elif args.mode == "reference":  # for reference, we are more interested in what is NOT Good
        query = etc.config.refEval(args.refGenome)
        good_ids = biotype_map.viewkeys() - sql_lib.get_query_ids(cur, query)  # actually not good
        out_good_bed_path, out_good_big_bed_path = get_bed_paths(args.outDir, "reference", args.refGenome)
        gp_dict = seq_lib.get_transcript_dict(args.annotationGp)
    elif args.mode == "transMap":
        good_ids = get_all_tm_good(cur, args.refGenome, args.genome)
        out_good_bed_path, out_good_big_bed_path = get_bed_paths(args.outDir, "transMap", args.genome)
        gp_dict = seq_lib.get_transcript_dict(args.targetGp)
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")
    with open(out_good_bed_path, "w") as outf:
        for aln_id, rec in gp_dict.iteritems():
            tx_id = psl_lib.strip_alignment_numbers(aln_id)
            if aln_id in good_ids:
                if biotype_map[tx_id] == "protein_coding":
                    bed = rec.get_bed(rgb=colors["coding"])
                    outf.write("".join(["\t".join(map(str, bed)), "\n"]))
                else:
                    bed = rec.get_bed(rgb=colors["noncoding"])
                    outf.write("".join(["\t".join(map(str, bed)), "\n"]))
            else:
                bed = rec.get_bed(rgb=colors["not_good"])
                outf.write("".join(["\t".join(map(str, bed)), "\n"]))
    make_big_bed(out_good_bed_path, args.sizes, out_good_big_bed_path)