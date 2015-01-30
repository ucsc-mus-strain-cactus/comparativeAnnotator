import os
import argparse
from itertools import izip_longest

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system
from lib.general_lib import FileType, DirType, FullPaths, classesInModule
import lib.sqlite_lib as sql_lib

import src.classifiers, src.details, src.attributes


#hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".2bit"
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
    parser.add_argument('--primaryKey', type=str, default="AlignmentId")
    return parser


def parseDir(genomes, targetDir, ext):
    pathDict = {}
    for g in genomes:
        path = os.path.join(targetDir, g + ext)
        if not os.path.exists(path):
            raise RuntimeError("{} does not exist".format(path))
        pathDict[g] = path
    return pathDict


def buildAnalyses(target, alnPslDict, seqTwoBitDict, refSeqTwoBit, geneCheckBedDict, 
            gencodeAttributeMap, genomes, annotationBed, outDir, primaryKeyColumn, refGenome):
    for genome in genomes:
        alnPsl = alnPslDict[genome]
        geneCheckBed = geneCheckBedDict[genome]
        seqTwoBit = seqTwoBitDict[genome]

        #find all user-defined classes in the classifiers module
        classifiers = classesInModule(src.classifiers)
        outClassify = os.path.join(outDir, "classify.db")
        initializeDb(outClassify, genome, classifiers, alnPsl, primaryKeyColumn)
        for c in classifiers:
            target.addChildTarget(c(genome, alnPsl, seqTwoBit, refSeqTwoBit,
                annotationBed, gencodeAttributeMap, geneCheckBed, outClassify, 
                refGenome, primaryKeyColumn))

        #find all user-defined classes in the details module
        details = classesInModule(src.details)        
        outDetails = os.path.join(outDir, "details.db")
        initializeDb(outDetails, genome, details, alnPsl, primaryKeyColumn)        
        for d in details:
            target.addChildTarget(d(genome, alnPsl, seqTwoBit, refSeqTwoBit,
                annotationBed, gencodeAttributeMap, geneCheckBed, outDetails, 
                refGenome, primaryKeyColumn))

        #find all user-defined classes in the attributes module     
        attributes = classesInModule(src.attributes)
        outAttributes = os.path.join(outDir, "attributes.db")
        initializeDb(outAttributes, genome, attributes, alnPsl, primaryKeyColumn)        
        for a in attributes:
            target.addChildTarget(a(genome, alnPsl, seqTwoBit, refSeqTwoBit,
                annotationBed, gencodeAttributeMap, geneCheckBed, outAttributes, 
                refGenome, primaryKeyColumn))


def initializeDb(dbPath, genome, classifiers, alnPsl, primaryKeyColumn):
    columnDefinitions = [[x.__name__, x._getType()] for x in classifiers]
    #find alignment IDs from PSLs (primary key for database)
    aIds = set(x.split()[9] for x in open(alnPsl))
    initializeSqlTable(dbPath, genome, columnDefinitions, primaryKeyColumn)
    initializeSqlRows(dbPath, genome, aIds, primaryKeyColumn)


def initializeSqlTable(db, genome, columns, primaryKey):
    with sql_lib.ExclusiveSqlConnection(db) as cur:
        sql_lib.initializeTable(cur, genome, columns, primaryKey)


def initializeSqlRows(db, genome, aIds, primaryKey):
    with sql_lib.ExclusiveSqlConnection(db) as cur:
        sql_lib.insertRows(cur, genome, primaryKey, [primaryKey], izip_longest(aIds, [None]))


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

    #find data files
    alnPslDict = parseDir(args.genomes, args.dataDir, alignment_ext)
    seqTwoBitDict = parseDir(args.genomes, args.dataDir, sequence_ext)
    geneCheckBedDict = parseDir(args.genomes, args.dataDir, gene_check_ext)

    refSeqTwoBit = os.path.join(args.dataDir, args.refGenome + ".2bit")
    if not os.path.exists(refSeqTwoBit):
        raise RuntimeError("Reference genome 2bit not present at {}".format(refSeqTwoBit))

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(alnPslDict, seqTwoBitDict, refSeqTwoBit, 
            geneCheckBedDict, args.gencodeAttributeMap, args.genomes, args.annotationBed, 
            args.outDir, args.primaryKey, args.refGenome))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from src.main import *
    main()