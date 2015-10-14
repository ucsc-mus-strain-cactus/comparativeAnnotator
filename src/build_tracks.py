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
            database(args.genome, db, db_path, tmp_dir)
    elif args.mode == "reference":
        for db in ["classify", "details"]:
            db_path = os.path.join(args.outDir, "{}.db".format(db))
            database(args.refGenome, db, db_path, tmp_dir, ref=True)
        attr_db_path = os.path.join(args.outDir, "attributes.db")
        ref_attr_table(args.refGenome, attr_db_path, args.gencodeAttributes)
    elif args.mode == "transMap":
        for db in ["classify", "details", "attributes"]:
            db_path = os.path.join(args.outDir, "{}.db".format(db))
            database(args.genome, db, db_path, tmp_dir)
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")
    target.setFollowOnTargetFn(build_tracks_wrapper, args=[args])


def database(genome, db, db_path, tmp_dir, ref=False):
    data_dict = {}
    mkdir_p(os.path.dirname(db_path))
    data_path = os.path.join(tmp_dir, db)
    for col in os.listdir(data_path):
        p = os.path.join(data_path, col)
        with open(p) as p_h:
            data_dict[col] = pickle.load(p_h)
    index_label = "AlignmentId" if ref else "TranscriptId"
    sql_lib.write_dict(data_dict, db_path, genome, index_label)


def ref_attr_table(ref_genome, db_path, attr_file):
    """
    This function is used to add an extra table in reference mode holding all of the basic attributes.
    Basically directly dumping the tsv into sqlite3.
    """
    sql_lib.write_csv(attr_file, db_path, ref_genome, index_col=3, sep="\t", index_label="TranscriptId")


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
    for f in classifier_tracks:
        query = f(genome)
        query_name = f.__name__
        target.addChildTargetFn(build_classifier_tracks, args=[query, query_name, genome, args])
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


def build_classifier_tracks(target, query, query_name, genome, args):
    con, cur = sql_lib.attach_databases(args.outDir, mode=args.mode)
    bed_recs = cur.execute(query)
    out_bed_path, out_big_bed_path = get_bed_paths(args.outDir, query_name, genome)
    with open(out_bed_path, "w") as outf:
        for recs in bed_recs:
            for rec in recs:
                if rec is not None:
                    outf.write(rec)
    make_big_bed(out_bed_path, args.sizes, out_big_bed_path)


def get_all_tm_good(cur, genome):
    """
    transMap Good varies depending on if the transcript is coding or noncoding. We will build a set of IDs for both.
    """
    biotypes = sql_lib.get_all_biotypes(cur, genome, gene_level=False)
    coding_query = etc.config.transMapEval(genome, biotype="protein_coding", good=True)
    coding_good = 
    non_coding_query = etc.config.transMapEval(genome, coding=False, good=True)
    coding_good_ids = sql_lib.get_query_ids(cur, coding_query)
    non_coding_good_ids = sql_lib.get_query_ids(cur, non_coding_query)
    coding_id_query = "SELECT AlignmentId FROM attributes.'{}' WHERE TranscriptType == 'protein_coding'".format(genome)
    coding_gencode_ids = {x[0] for x in cur.execute(coding_id_query).fetchall()}
    coding_ok = coding_gencode_ids & coding_good_ids
    non_coding_id_query = "SELECT AlignmentId FROM attributes.'{}' WHERE TranscriptType != 'protein_coding'".format(genome)
    non_coding_gencode_ids = {x[0] for x in cur.execute(non_coding_id_query).fetchall()}
    non_coding_ok = non_coding_gencode_ids & non_coding_good_ids
    return non_coding_ok, coding_ok


def build_good_track(target, args):
    """
    Builds a specific track of Good transcripts for the current mode.
    """
    colors = {"coding": "100,209,61", "noncoding": "209,61,115"}
    con, cur = sql_lib.attach_databases(args.outDir, args.mode)
    if args.mode == "augustus":
        query = etc.config.augustusEval(args.genome)
        good_ids = sql_lib.get_query_ids(cur, query)
        out_bed_path, out_big_bed_path = get_bed_paths(args.outDir, "augustusOk", args.genome)
        gp_dict = seq_lib.get_transcript_dict(args.augustusGp)
    elif args.mode == "reference":
        query = etc.config.refEval(args.refGenome)
        good_ids = sql_lib.get_query_ids(cur, query)
        out_bed_path, out_big_bed_path = get_bed_paths(args.outDir, "referenceOk", args.refGenome)
        gp_dict = seq_lib.get_transcript_dict(args.annotationGp)
    else:
        non_coding_ok, coding_ok = get_all_tm_good(cur, args.genome)
        good_ids = non_coding_ok & coding_ok
        out_bed_path, out_big_bed_path = get_bed_paths(args.outDir, "transMapOk", args.genome)
        gp_dict = seq_lib.get_transcript_dict(args.targetGp)
    with open(out_bed_path, "w") as outf:
        for aln_id, rec in gp_dict.iteritems():
            if aln_id in good_ids:
                if args.mode == "transMap" and aln_id in non_coding_ok:
                    bed = rec.get_bed(rgb=colors["noncoding"])
                    outf.write("".join(["\t".join(map(str, bed)), "\n"])) 
                else:
                    bed = rec.get_bed(rgb=colors["coding"])
                    outf.write("".join(["\t".join(map(str, bed)), "\n"]))
    make_big_bed(out_bed_path, args.sizes, out_big_bed_path)