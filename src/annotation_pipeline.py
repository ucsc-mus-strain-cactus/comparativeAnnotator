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
import src.alignment_classifiers
import src.augustus_classifiers
import src.attributes
import etc.config

from src.build_tracks import database_wrapper

__author__ = "Ian Fiddes"


def parse_args():
    """
    Builds an argument parser for this run
    """
    parent_parser = argparse.ArgumentParser()
    subparsers = parent_parser.add_subparsers(title="Modes", description="Execution Modes", dest="mode")
    tm_parser = subparsers.add_parser('transMap')
    ref_parser = subparsers.add_parser('reference')
    aug_parser = subparsers.add_parser('augustus')
    # common arguments
    for parser in [ref_parser, aug_parser, tm_parser]:
        parser.add_argument('--outDir', type=str, required=True)
        parser.add_argument('--refGenome', type=str, required=True)
        parser.add_argument('--refFasta', type=str, required=True)
        parser.add_argument('--sizes', required=True)
        parser.add_argument('--annotationGp', required=True)
    # transMap specific options
    for parser in [aug_parser, tm_parser]:
        parser.add_argument('--gencodeAttributes', required=True)
        parser.add_argument('--genome', required=True)
        parser.add_argument('--psl', required=True)
        parser.add_argument('--targetGp', required=True)
        parser.add_argument('--fasta', required=True)
    # Augustus specific options
    aug_parser.add_argument('--augustusGp', required=True)
    Stack.addJobTreeOptions(parent_parser)  # add jobTree options
    args = parent_parser.parse_args()
    return args


def run_ref_classifiers(args, target, tmp_dir):
    ref_classifiers = classes_in_module(src.classifiers)
    for classifier in ref_classifiers:
        target.addChildTarget(classifier(args.refFasta, args.annotationGp, args.refGenome, tmp_dir))


def run_tm_classifiers(args, target, tmp_dir):
    tm_classifiers = classes_in_module(src.alignment_classifiers)
    for classifier in tm_classifiers:
        target.addChildTarget(classifier(args.refFasta, args.annotationGp, args.refGenome, tmp_dir, args.genome,
                                         args.psl, args.fasta, args.targetGp))


def run_aug_classifiers(args, target, tmp_dir):
    aug_classifiers = classes_in_module(src.augustus_classifiers)
    for classifier in aug_classifiers:
        target.addChildTarget(classifier(args.refFasta, args.annotationGp, args.refGenome, tmp_dir, args.genome,
                                         args.psl, args.fasta, args.targetGp, args.augustusGp))


def build_analyses(target, args):
    """
    Wrapper function that will call all classifiers. Each classifier will dump its results to disk as a pickled dict.
    Calls database_wrapper to load these into a sqlite3 database.
    """
    tmp_dir = target.getGlobalTempDir()
    if args.mode in ["reference", "transMap", "augustus"]:
        run_ref_classifiers(args, target, tmp_dir)
    if args.mode in ["transMap, augustus"]:
        run_tm_classifiers(args, target, tmp_dir)
    if args.mode == "augustus":
        run_aug_classifiers(args, target, tmp_dir)
    # merge the resulting pickled files into sqlite databases and construct BED tracks
    target.setFollowOnTargetFn(database_wrapper, memory=8 * (1024 ** 3), args=[args, tmp_dir])


def main():
    args = parse_args()
    i = Stack(Target.makeTargetFn(build_analyses, memory=8 * (1024 ** 3), args=[args])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotation_pipeline import *
    main()
