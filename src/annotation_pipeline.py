"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""

import argparse
import os
import subprocess
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
from lib.general_lib import classes_in_module, mkdir_p

import src.classifiers
import src.augustus_classifiers
import src.attributes
import etc.config

from src.build_tracks import database_wrapper

__author__ = "Ian Fiddes"


def parse_args():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str, required=True)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--annotationGp', required=True)
    parser.add_argument('--psl', required=True)
    parser.add_argument('--gp', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--refFasta', type=str, required=True)
    parser.add_argument('--sizes', required=True)
    parser.add_argument('--gencodeAttributes', required=True)
    parser.add_argument('--outDir', type=str, required=True)
    parser.add_argument('--augustus', action="store_true")
    parser.add_argument('--augustusGp')
    Stack.addJobTreeOptions(parser)  # add jobTree options
    args = parser.parse_args()
    if args.augustus and args.augustusGp is None:
        raise RuntimeError("Error: augustus mode activated but no augustusGp provided.")
    elif args.augustus is False and args.augustusGp is not None:
        raise RuntimeError("Error: augustusGp provided but augustus flag not set.")
    return args


def build_analyses(target, ref_genome, genome, annotation_gp, psl, gp, fasta, ref_fasta, sizes, gencode_attributes,
                   out_dir, augustus, augustus_gp):
    """
    Wrapper function that will call all classifiers. Each classifier will dump its results to disk as a pickled dict.
    Calls database_wrapper to load these into a sqlite3 database.
    """
    tmp_dir = target.getGlobalTempDir()
    if augustus is True:
        # find all user-defined classes in the categories of analyses
        augustus_classifiers = classes_in_module(src.augustus_classifiers)
        for classifier in augustus_classifiers:
            target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                             ref_genome, tmp_dir, augustus_gp))
        target.setFollowOnTargetFn(database_wrapper, memory=8 * (1024 ** 3),
                                   args=[out_dir, genome, sizes, augustus_gp, augustus, tmp_dir])
    else:
        # find all user-defined classes in the categories of analyses
        classifiers = classes_in_module(src.classifiers) + classes_in_module(src.attributes)
        for classifier in classifiers:
            target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                             ref_genome, tmp_dir))
        # merge the resulting pickled files into sqlite databases and construct BED tracks
        target.setFollowOnTargetFn(database_wrapper, memory=8 * (1024 ** 3),
                                   args=[out_dir, genome, sizes, gp, augustus, tmp_dir])


def main():
    args = parse_args()
    i = Stack(Target.makeTargetFn(build_analyses, memory=8 * (1024 ** 3),
                                  args=[args.refGenome, args.genome, args.annotationGp, args.psl, args.gp, args.fasta, 
                                        args.refFasta, args.sizes, args.gencodeAttributes, args.outDir, args.augustus,
                                        args.augustusGp])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotation_pipeline import *
    main()
