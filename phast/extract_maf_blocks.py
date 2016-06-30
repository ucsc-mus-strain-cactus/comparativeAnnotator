"""
jobTree wrapper for extracting MAF blocks. Extracts each block in the reference that is conserved.
Meant to be used downstream of the main phastCons pipeline in the process of detecting accelerated regions.
"""

import os
import argparse
import itertools
import shutil
import random
import cPickle as pickle
from ete3 import Tree
from collections import OrderedDict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, TempFileTree
from phast.phast_functions import *
from phast.phast_subset import subset_hal_pipeline


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hal', help='HAL alignment file.')
    parser.add_argument('ref_genome', help='Reference genome.')
    parser.add_argument('conserved_bed', help='Conserved BED file.')
    parser.add_argument('--target_genomes', nargs='+', required=True, help='target genomes')
    parser.add_argument('--mafDir', help='output base dir.', default='mafBlocks')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def extract_maf_wrapper(target, args):
    """
    Main pipeline wrapper. Calls out to hal2maf once for each region in args.conserved_bed
    """
    bed_recs = [x.split()[:3] for x in open(args.conserved_bed)]
    os.makedirs(args.mafDir)
    file_tree = TempFileTree(args.mafDir)
    for chrom, start, stop in bed_recs:
        p = file_tree.getTempFile(suffix='_' + '_'.join([chrom, start, stop]))
        target.addChildTargetFn(hal2maf, args=(args, chrom, start, stop, p))


def hal2maf(target, args, chrom, start, stop, p):
    """
    Runs hal2maf on the region defined by the rec.
    """
    size = int(stop) - int(start)
    cmd = 'hal2maf --targetGenomes {} --unique --noDupes --noAncestors --refGenome {} --refSequence {} --start {} --length {} {} {}'
    cmd = cmd.format(','.join(args.target_genomes), args.ref_genome, chrom, start, size, args.hal, p)
    system(cmd)


if __name__ == '__main__':
    from phast.extract_maf_blocks import *
    args = parse_args()
    s = Stack(Target.makeTargetFn(extract_maf_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
