"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""
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
        parser.add_argument('--outDir', type=str, required=True)
        parser.add_argument('--refGenome', type=str, required=True)
        parser.add_argument('--refFasta', type=str, required=True)
        parser.add_argument('--sizes', required=True)
        parser.add_argument('--annotationGp', required=True)
        parser.add_argument('--gencodeAttributes', required=True)
        Stack.addJobTreeOptions(parser)  # add jobTree options
    # transMap specific options
    for parser in [aug_parser, tm_parser]:
        parser.add_argument('--genome', required=True)
        parser.add_argument('--psl', required=True)
        parser.add_argument('--refPsl', required=True)
        parser.add_argument('--targetGp', required=True)
        parser.add_argument('--fasta', required=True)
    # Augustus specific options
    aug_parser.add_argument('--augustusGp', required=True)
    args = parent_parser.parse_args()
    assert args.mode in ["transMap", "reference", "augustus"]
    return args


def build_analyses(target, args):
    """
    Wrapper function that will call all classifiers. Each classifier will dump its results to disk as a pickled dict.
    Calls database_wrapper to load these into a sqlite3 database.
    """
    tmp_dir = target.getGlobalTempDir()
    if args.mode == "reference":
        run_ref_classifiers(args, target, tmp_dir)
    elif args.mode == "transMap":
        run_tm_classifiers(args, target, tmp_dir)
    elif args.mode == "augustus":
        run_aug_classifiers(args, target, tmp_dir)
    else:
        raise RuntimeError("Somehow your argparse object does not contain a valid mode.")


def main():
    args = parse_args()
    i = Stack(Target.makeTargetFn(build_analyses, memory=8 * (1024 ** 3), args=[args])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from comparativeAnnotator.annotation_pipeline import *
    main()
