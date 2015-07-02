import os
import argparse
from itertools import izip

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger
from lib.general_lib import classesInModule
import lib.sqlite_lib as sql_lib

import src.classifiers
import src.augustusClassifiers
import src.attributes
from src.constructDatabases import ConstructDatabases, ConstructAugustusDatabases
from src.buildTracks import BuildTracks


def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str, required=True)
    parser.add_argument('--genomes', nargs="+", required=True)
    parser.add_argument('--annotationGp', required=True)
    parser.add_argument('--psls', nargs='+', required=True)
    parser.add_argument('--gps', nargs='+', required=True)
    parser.add_argument('--augustusGps', nargs='+', required=True)
    parser.add_argument('--fastas', nargs='+', required=True)
    parser.add_argument('--refFasta', type=str, required=True)
    parser.add_argument('--sizes', nargs='+', required=True)
    parser.add_argument('--gencodeAttributeMap', required=True)
    parser.add_argument('--outDir', type=str, required=True)
    parser.add_argument('--primaryKeyColumn', type=str, default="AlignmentId")
    return parser


def buildAnalyses(target, psls, fastas, refFasta, gps, gencodeAttributeMap, genomes,
                  annotationGp, outDir, refGenome, primaryKeyColumn, sizes, augGps):
    # find all user-defined classes in the categories of analyses
    classifiers = classesInModule(src.classifiers)
    augustusClassifiers = classesInModule(src.augustusClassifiers)
    attributes = classesInModule(src.attributes)
    for genome, psl, bed, fasta, augGp in izip(genomes, psls, gps, fastas, augGps):
        for c in classifiers:
            target.addChildTarget(c(genome, psl, fasta, refFasta, annotationGp, gencodeAttributeMap,
                                    bed, refGenome, primaryKeyColumn, target.getGlobalTempDir()))
        for a in attributes:
            target.addChildTarget(a(genome, psl, fasta, refFasta, annotationGp, gencodeAttributeMap,
                                    bed, refGenome, primaryKeyColumn, target.getGlobalTempDir()))
        for aug in augustusClassifiers:
            target.addChildTarget(aug(genome, psl, fasta, refFasta, annotationGp, gencodeAttributeMap,
                                      bed, refGenome, primaryKeyColumn, target.getGlobalTempDir(), augGp))
        # merge the resulting pickled files into sqlite databases
    target.setFollowOnTargetFn(databaseWrapper, args=(outDir, genomes, psls, primaryKeyColumn, sizes, gps, augGps,
                                                      annotationGp))


def databaseWrapper(target, outDir, genomes, psls, primaryKeyColumn, sizes, gps, augGps, annotationGp):
    target.addChildTarget(ConstructDatabases(outDir, target.getGlobalTempDir(), genomes, psls, primaryKeyColumn))
    target.addChildTarget(ConstructAugustusDatabases(outDir, target.getGlobalTempDir(), genomes, augGps,
                                                     primaryKeyColumn))
    target.setFollowOnTarget(BuildTracks(outDir, genomes, primaryKeyColumn, sizes, gps, annotationGp))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    assert len(args.psls) == len(args.fastas) == len(args.genomes) == len(args.gps) == len(args.sizes)
    for genome, psl, gp, fasta, size in izip(args.genomes, args.psls, args.gps, args.fastas, args.sizes):
        assert all([genome in x for x in [psl, gp, fasta, size]])

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)

    if not os.path.exists(args.refFasta):
        raise RuntimeError("Reference genome fasta not present at {}".format(args.refFasta))
    elif not os.path.exists(args.annotationGp):
        raise RuntimeError("Annotation bed not present at {}".format(args.annotationGp))

    for x in args.psls:
            if not os.path.exists(x):
                raise RuntimeError("PSL not present at {}".format(x))

    for x in args.gps:
            if not os.path.exists(x):
                raise RuntimeError("BED not present at {}".format(x))

    for x in args.fastas:
            if not os.path.exists(x):
                raise RuntimeError("Fasta not present at {}".format(x))

    for x in args.sizes:
            if not os.path.exists(x):
                raise RuntimeError("chrom.sizes not present at {}".format(x))

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(sorted(args.psls), sorted(args.fastas), args.refFasta,
                                                       sorted(args.gps), args.gencodeAttributeMap, sorted(args.genomes),
                                                       args.annotationGp, args.outDir, args.refGenome, 
                                                       args.primaryKeyColumn, sorted(args.sizes), args.augustusGps,
                                                       ))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.annotationPipeline import *
    main()
