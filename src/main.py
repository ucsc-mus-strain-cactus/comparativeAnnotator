import os
import argparse

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system
from lib.general_lib import FileType, DirType, FullPaths

from src.classifierMaster import ClassifierMaster

#hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".2bit"
gene_check_ext = ".bed"

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
    parser.add_argument('--outDir', type=str, default="output/")
    parser.add_argument('--primaryKey', type=str, default="AlignmentID")
    parser.add_argument('--overwriteDb', action="store_true")
    parser.add_argument('--mergedDb', type=str, default="results.db")
    return parser


def parse_dir(genomes, targetDir, ext):
    pathDict = {}
    for g in genomes:
        path = os.path.join(targetDir, g + ext)
        if not os.path.exists(path):
            raise RuntimeError("{} does not exist".format(path))
        pathDict[g] = path
    return pathDict


def build_analysis(target, alnPslDict, seqTwoBitDict, refSeqTwoBit, geneCheckBedDict, 
            gencodeAttributeMap, genomes, annotationBed, outDir, primaryKeyColumn, refGenome):
    for genome in genomes:
        alnPsl = alnPslDict[genome]
        geneCheckBed = geneCheckBedDict[genome]
        seqTwoBit = seqTwoBitDict[genome]
        target.addChildTarget(ClassifierMaster(genome, alnPsl, seqTwoBit, refSeqTwoBit, annotationBed,
                gencodeAttributeMap, geneCheckBed, outDir, refGenome, primaryKeyColumn))


def initialize_sql_columns(genome, outDir, primaryKeyColumn):
    outDb = os.path.join(outDir, genome + ".db")
    con = sql.connect(outDb)
    columns = [[x.__name__, x.__type__()] for x in classifiers]
    with con:
        initializeTable(con.cursor(), genome, columns, primaryKeyColumn)


def initialize_sql_rows(genome, outDir, alnPsl, primaryKeyColumn):
    outDb = os.path.join(outDir, genome + ".db")
    con = sql.connect(outDb)
    alnIds = set(x.split()[9] for x in open(alnPsl))
    for alnId in alnIds:
        insertRow(con.cursor(), genome, primaryKeyColumn, alnId)


def merge_databases(outDir, mergedDb, genomes):
    dbs = [os.path.join(outDir, x + ".db") for x in genomes]
    for db in dbs:
        system("sqlite3 {} .dump | sqlite3 {}".format(db, mergedDb))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)

    if args.overwriteDb is True:
        if os.path.exists(args.mergedDb):
            os.remove(args.mergedDb)
        for g in args.genomes:
            if os.path.exists(os.path.join(args.outDir, g + ".db")):
                os.remove(os.path.join(args.outDir, g + ".db"))

    alnPslDict = parse_dir(args.genomes, args.dataDir, alignment_ext)
    seqTwoBitDict = parse_dir(args.genomes, args.dataDir, sequence_ext)
    geneCheckBedDict = parse_dir(args.genomes, args.dataDir, gene_check_ext)

    refSeqTwoBit = os.path.join(args.dataDir, args.refGenome + ".2bit")
    if not os.path.exists(refSeqTwoBit):
        raise RuntimeError("Reference genome 2bit not present at {}".format(refSeqTwoBit))

    i = Stack(Target.makeTargetFn(build_analysis, args=(alnPslDict, seqTwoBitDict, refSeqTwoBit, 
            geneCheckBedDict, args.gencodeAttributeMap, args.genomes, args.annotationBed, args.outDir, 
            args.primaryKey, args.refGenome))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

    merge_databases(args.outDir, args.mergedDb, args.genomes)


if __name__ == '__main__':
    from src.main import *
    main()