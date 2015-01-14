import os
import argparse
import sqlite3 as sql
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions, system
from lib.sqlite_lib import initializeTable, insertRow
from lib.general_lib import FileType, DirType, FullPaths

#basic_attributes has many basic classes stuck together that are not really classifiers
from src.basic_attributes import *
#psl_attributes takes attributes from the psl and makes columns, mainly genomic positions
from src.psl_attributes import *

from src.unknown_bases import UnknownBases
from src.end_stop import EndStop
from src.begin_start import BeginStart
from src.in_frame_stop import InFrameStop
from src.bad_frame import BadFrame
from src.no_cds import NoCds
from src.cds_gap import CdsGap
from src.cds_mult_3_gap import CdsMult3Gap
from src.utr_gap import UtrGap
from src.cds_unknown_splice import CdsUnknownSplice
from src.cds_non_canon_splice import CdsNonCanonSplice
from src.utr_unknown_splice import UtrUnknownSplice
from src.utr_non_canon_splice import UtrNonCanonSplice
from src.number_scaffold_gap import NumberScaffoldGap
from src.alignment_identity import AlignmentIdentity
from src.alignment_coverage import AlignmentCoverage
from src.alignment_partial_map import AlignmentPartialMap
from src.alignment_abuts_right import AlignmentAbutsRight
from src.alignment_abuts_left import AlignmentAbutsLeft

#classifiers we are currently working with
classifiers = [EndStop, UnknownBases, BeginStart, InFrameStop, BadFrame, NoCds, CdsMult3Gap, 
        UtrGap, CdsUnknownSplice, CdsNonCanonSplice, UtrUnknownSplice, UtrNonCanonSplice,
        NumberScaffoldGap, AlignmentIdentity, AlignmentCoverage, AlignmentPartialMap, 
        AlignmentAbutsRight, AlignmentAbutsLeft]

#add in all of the basic attribute columns
classifiers = classifiers + [TranscriptID, GeneID, GeneName, GeneType, TranscriptType]
#add in all of the psl attribute columns
classifiers = classifiers + [SourceChrom, SourceStart, SourceStop, SourceStrand,
                            DestChrom, DestStart, DestStop, DestStrand]


#hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".2bit"
gene_check_ext = ".bed"
#gene_check_details_ext = ".coding-gene-check-details.bed"

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


def build_analysis(target, alnPslDict, seqTwoBitDict, geneCheckBedDict, gencodeAttributeMap,
            genomes, annotationBed, outDir, primaryKeyColumn, refGenome):
    for genome in genomes:
        alnPsl, seqFasta = alnPslDict[genome], seqTwoBitDict[genome]
        geneCheckBed = geneCheckBedDict[genome]
        initialize_sql_columns(genome, outDir, primaryKeyColumn)
        for classifier in classifiers:
            target.addChildTarget(classifier(genome, alnPsl, seqFasta, annotationBed,
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

    logger.info("Building paths to the required files")
    alnPslDict = parse_dir(args.genomes, args.dataDir, alignment_ext)
    seqTwoBitDict = parse_dir(args.genomes, args.dataDir, sequence_ext)
    geneCheckBedDict = parse_dir(args.genomes, args.dataDir, gene_check_ext)
    #geneCheckBedDetailsDict = parse_dir(args.genomes, args.geneCheckDir, gene_check_details_ext)

    refSequence = os.path.join(args.dataDir, args.refGenome + ".2bit")
    if not os.path.exists(refSequence):
        raise RuntimeError("Reference genome 2bit not present at {}".format(refSequence))
    args.refSequence = refSequence

    i = Stack(Target.makeTargetFn(build_analysis, args=(alnPslDict, seqTwoBitDict, geneCheckBedDict, 
            args.gencodeAttributeMap, args.genomes, args.annotationBed, args.outDir, args.primaryKey, 
            args.refGenome))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

    merge_databases(args.outDir, args.mergedDb, args.genomes)


if __name__ == '__main__':
    from src.main import *
    main()