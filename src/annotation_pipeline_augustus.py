"""
This is the main driver script for comparativeAnnotator in AugustusTMR mode.
"""

import argparse

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import TempFileTree

from lib.general_lib import classes_in_module

import src.augustus_classifiers
from src.construct_databases import ConstructDatabases
from src.build_tracks import BuildTracks

__author__ = "Ian Fiddes"


def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str, required=True)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--annotationGp', required=True)
    parser.add_argument('--psl', required=True)
    parser.add_argument('--gp', required=True)
    parser.add_argument('--augustusGp', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--refFasta', type=str, required=True)
    parser.add_argument('--sizes', required=True)
    parser.add_argument('--gencodeAttributeMap', required=True)
    parser.add_argument('--outDir', type=str, required=True)
    return parser


def build_analyses(target, ref_genome, genome, annotation_gp, psl, gp, aug_gp, fasta, ref_fasta, sizes,
                   gencode_attributes, out_dir):
    # find all user-defined classes in the categories of analyses
    out_file_tree = TempFileTree(target.getGlobalTempDir())
    classifiers = classes_in_module(src.augustus_classifiers)
    for classifier in classifiers:
        target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                         ref_genome, out_file_tree, aug_gp))
        # merge the resulting pickled files into sqlite databases and construct BED tracks
    target.setFollowOnTargetFn(database, memory=8 * (1024 ** 3),
                               args=(out_dir, genome, psl, sizes, gp, annotation_gp, out_file_tree))


def database(target, out_dir, genome, psl, sizes, gp, annotation_gp, out_file_tree):
    target.addChildTarget(ConstructDatabases(out_dir, out_file_tree, genome, psl))
    target.setFollowOnTarget(BuildTracks(out_dir, genome, sizes, gp, annotation_gp))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    i = Stack(Target.makeTargetFn(build_analyses, args=(args.refGenome, args.genome, args.annotationGp, args.psl,
                                                        args.gp, args.augustusGp, args.fasta, args.refFasta, args.sizes,
                                                        args.gencodeAttributes, args.outDir))).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotation_pipeline_augustus import *
    main()
