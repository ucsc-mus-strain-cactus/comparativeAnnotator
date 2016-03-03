"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""
import sys
import argparse
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from comparativeAnnotator.classify_driver import run_ref_classifiers, run_tm_classifiers, run_aug_classifiers

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
        parser.add_argument('--db', type=str, required=True)
        parser.add_argument('--ref_genome', type=str, required=True)
        parser.add_argument('--ref_fasta', type=str, required=True)
        parser.add_argument('--sizes', required=True)
        parser.add_argument('--annotation_gp', required=True)
        parser.add_argument('--gencode_attributes', required=True)
        Stack.addJobTreeOptions(parser)  # add jobTree options
    # transMap specific options
    for parser in [aug_parser, tm_parser]:
        parser.add_argument('--genome', required=True)
        parser.add_argument('--psl', required=True)
        parser.add_argument('--ref_psl', required=True)
        parser.add_argument('--target_gp', required=True)
        parser.add_argument('--fasta', required=True)
    # Augustus specific options
    aug_parser.add_argument('--augustus_gp', required=True)
    args = parent_parser.parse_args()
    assert args.mode in ["transMap", "reference", "augustus"]
    return args


def comp_ann_driver(target, args):
    """
    Wrapper function that will call all classifiers. Each classifier will dump its results to disk as a pickled dict.
    Calls database_wrapper to load these into a sqlite3 database.
    """
    if args.mode == "reference":
        run_ref_classifiers(target, args)
    elif args.mode == "transMap":
        run_tm_classifiers(target, args)
    elif args.mode == "augustus":
        run_aug_classifiers(target, args)
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")


def main(args):
    i = Stack(Target.makeTargetFn(comp_ann_driver, args=[args])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from comparativeAnnotator.annotation_pipeline import *
    args = parse_args()
    main(args)
