import os
import argparse

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger
from lib.general_lib import FileType, DirType, FullPaths, classesInModule
import lib.sqlite_lib as sql_lib

import src.classifiers, src.details, src.attributes
from src.constructDatabases import ConstructDatabases
from src.buildTracks import BuildTracks


# hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".fa"
gene_check_ext = ".gene-check.bed"


def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str)
    parser.add_argument('--genomes', nargs="+")
    parser.add_argument('--annotationBed', type=FileType)
    parser.add_argument('--dataDir', type=DirType, action=FullPaths)
    parser.add_argument('--gencodeAttributeMap', type=FileType)
    parser.add_argument('--outDir', type=str, default="./output/", action=FullPaths)
    parser.add_argument('--primaryKeyColumn', type=str, default="AlignmentId")
    return parser


def parseDir(genomes, targetDir, ext):
    pathDict = {}
    for g in genomes:
        path = os.path.join(targetDir, g + ext)
        if not os.path.exists(path):
            raise RuntimeError("{} does not exist".format(path))
        pathDict[g] = path
    return pathDict


def buildAnalyses(target, alnPslDict, fastaDict, refSeqTwoBit, geneCheckBedDict, gencodeAttributeMap, genomes,
                  annotationBed, outDir, refGenome, primaryKeyColumn, dataDir):
    # find all user-defined classes in the three categories of analyses
    classifiers = classesInModule(src.classifiers)
    details = classesInModule(src.details)
    attributes = classesInModule(src.attributes)
    for genome in genomes:
        alnPsl = alnPslDict[genome]
        geneCheckBed = geneCheckBedDict[genome]
        fasta = fastaDict[genome]
        # set child targets for every classifier-genome pair
        for c in classifiers:
            target.addChildTarget(
                c(genome, alnPsl, fasta, refSeqTwoBit, annotationBed, gencodeAttributeMap, geneCheckBed, refGenome,
                  primaryKeyColumn, outDir, "classify"))
        for d in details:
            target.addChildTarget(
                d(genome, alnPsl, fasta, refSeqTwoBit, annotationBed, gencodeAttributeMap, geneCheckBed, refGenome,
                  primaryKeyColumn, outDir, "details"))
        for a in attributes:
            target.addChildTarget(
                a(genome, alnPsl, fasta, refSeqTwoBit, annotationBed, gencodeAttributeMap, geneCheckBed, refGenome,
                  primaryKeyColumn, outDir, "attributes"))
            #merge the resulting pickled files into sqlite databases
    target.setFollowOnTargetFn(databaseWrapper, args=(
    outDir, genomes, classifiers, details, attributes, alnPslDict, primaryKeyColumn, 
    dataDir, geneCheckBedDict, annotationBed))


def databaseWrapper(target, outDir, genomes, classifiers, details, attributes, alnPslDict, primaryKeyColumn,
                     dataDir, geneCheckBedDict, annotationBed):
    target.addChildTarget(
        ConstructDatabases(outDir, genomes, classifiers, details, attributes, alnPslDict, primaryKeyColumn))
    target.setFollowOnTarget(BuildTracks(outDir, genomes, classifiers, details, attributes, primaryKeyColumn,
                      dataDir, geneCheckBedDict, annotationBed))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)

    for x in os.listdir(args.outDir):
        if x.endswith(".db"):
            os.remove(os.path.join(args.outDir, x))

    # find data files
    alnPslDict = parseDir(args.genomes, args.dataDir, alignment_ext)
    fastaDict = parseDir(args.genomes, args.dataDir, sequence_ext)
    geneCheckBedDict = parseDir(args.genomes, args.dataDir, gene_check_ext)

    refSeqTwoBit = os.path.join(args.dataDir, args.refGenome + ".2bit")
    if not os.path.exists(refSeqTwoBit):
        raise RuntimeError("Reference genome 2bit not present at {}".format(refSeqTwoBit))

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(alnPslDict, fastaDict, refSeqTwoBit,
                                                       geneCheckBedDict, args.gencodeAttributeMap, args.genomes,
                                                       args.annotationBed,
                                                       args.outDir, args.refGenome, args.primaryKeyColumn,
                                                       args.dataDir))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.main import *

    main()
