"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""

import argparse

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from lib.general_lib import classes_in_module

import src.classifiers
import src.augustus_classifiers
import src.attributes
from src.construct_databases import ConstructDatabases
from src.build_tracks import BuildTracks

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
    return args


def build_analyses(target, ref_genome, genome, annotation_gp, psl, gp, fasta, ref_fasta, sizes, gencode_attributes,
                   out_dir, augustus, augustus_gp):
    # find all user-defined classes in the categories of analyses
    if augustus is True:
        augustus_classifiers = classes_in_module(src.augustus_classifiers)
        for classifier in augustus_classifiers:
            target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                             ref_genome, augustus_gp))
    else:
        classifiers = classes_in_module(src.classifiers) + classes_in_module(src.attributes)
        for classifier in classifiers:
            target.addChildTarget(classifier(genome, psl, fasta, ref_fasta, annotation_gp, gencode_attributes, gp,
                                             ref_genome))
        # merge the resulting pickled files into sqlite databases and construct BED tracks
    target.setFollowOnTargetFn(database, memory=8 * (1024 ** 3),
                               args=[out_dir, genome, psl, sizes, gp, annotation_gp])


def database(target, out_dir, genome, psl, sizes, gp, annotation_gp):
    target.addChildTarget(ConstructDatabases(out_dir, genome, psl))
    target.setFollowOnTarget(BuildTracks(out_dir, genome, sizes, gp, annotation_gp))


def main():
    args = parse_args()
    i = Stack(Target.makeTargetFn(build_analyses, args=[args.refGenome, args.genome, args.annotationGp, args.psl,
                                                        args.gp, args.fasta, args.refFasta, args.sizes,
                                                        args.gencodeAttributes, args.outDir, args.augustus,
                                                        args.augustutsGp])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotation_pipeline import *
    main()
