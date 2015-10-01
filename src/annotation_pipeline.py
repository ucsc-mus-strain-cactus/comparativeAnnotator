import os
import argparse
import itertools

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import TempFileTree

from lib.general_lib import classes_in_module
import lib.sqlite_lib as sql_lib

import src.classifiers
import src.attributes
from src.constructDatabases import ConstructDatabases
from src.buildTracks import BuildTracks


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
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--refFasta', type=str, required=True)
    parser.add_argument('--sizes', required=True)
    parser.add_argument('--gencodeAttributes', required=True)
    parser.add_argument('--outDir', type=str, required=True)
    return parser


def build_analyses(target, ref_genome, genome, annotation_gp, psl, gp, fasta, ref_fasta, sizes, gencode_attributes,
                   out_dir):
    # find all user-defined classes in the categories of analyses
    out_file_tree = TempFileTree(target.getGlobalTempDir())
    classifiers = classes_in_module(src.classifiers)
    attributes = classes_in_module(src.attributes)
    for classifier in classifiers:
        target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                         ref_genome, out_file_tree))
    for attribute in attributes:
        target.addChildTarget(attribute(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                        ref_genome, out_file_tree))
        # merge the resulting pickled files into sqlite databases and construct BED tracks
    target.setFollowOnTargetFn(database, memory=8 * (1024 ** 3),
                               args=(out_dir, genome, psl, sizes, gp, annotation_gp, out_file_tree))


def database(target, out_dir, genome, psl, sizes, gp, annotation_gp):
    target.addChildTarget(ConstructDatabases(outDir, target.getGlobalTempDir(), genomes, psls))
    target.setFollowOnTarget(BuildTracks(outDir, genomes, sizes, gps, annotationGp))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()

    i = Stack(Target.makeTargetFn(build_analyses, args=(args.refGenome, args.genome, args.annotationGp, args.psl,
                                                       args.gp, args.fasta, args.refFasta, args.sizes,
                                                       args.gencodeAttributes, args.outDir))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotationPipeline import *
    main()
