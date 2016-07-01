"""
Run acceleration tests on each region in the conserved_bed.
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
    parser.add_argument('--hal', help='HAL alignment file.', required=True)
    parser.add_argument('--ref_genome', help='Reference genome.', required=True)
    parser.add_argument('--conserved_bed', help='Conserved BED file.', required=True)
    parser.add_argument('--conserved_model', help='Original conserved model', required=True)
    parser.add_argument('--non_conserved_model', help='Original nonconserved model.', required=True)
    parser.add_argument('--out_bed', help='output BED', required=True)
    parser.add_argument('--target_genomes', nargs='+', required=True, help='target genomes')
    parser.add_argument('--accelerated_genomes', nargs='+', required=True, help='target genomes')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def extract_maf_wrapper(target, args):
    """
    Main pipeline wrapper. Calls out to hal2maf once for each region in args.conserved_bed
    """
    bed_recs = [x.split()[:3] for x in open(args.conserved_bed)]
    result_dir = target.getGlobalTempDir()
    result_tree = TempFileTree(result_dir)
    for chrom, start, stop in bed_recs:
        result_path = result_tree.getTempFile(suffix='_' + '_'.join([chrom, start, stop]))
        target.addChildTargetFn(extract_and_calculate, args=(args, chrom, start, stop, result_path))
    target.setFollowOnTargetFn(cat_results, args=(args, result_tree.listFiles()))


def rename_model(model, accelerated_genomes):
    """rename the ancestor of accelerated_genomes. Necessary to do 2nd fitting"""
    lines = open(model).readlines()
    t = Tree(lines[-1].split('TREE: ')[1], format=1)
    anc = t.get_common_ancestor(accelerated_genomes)
    anc.name = 'Anc'
    with open(model, 'w') as outf:
        for l in lines[:-1]:
            outf.write(l)
        outf.write('TREE: ' + t.write(format=1) + '\n')


def extract_and_calculate(target, args, chrom, start, stop, result_path):
    """
    Runs the fitting.
    """
    # extract region MAF
    maf_path = os.path.join(target.getLocalTempDir(), 'region.maf')
    size = int(stop) - int(start)
    cmd = 'hal2maf --noAncestors --targetGenomes {} --refGenome {} --refSequence {} --start {} --length {} {} {}'
    cmd = cmd.format(','.join(args.target_genomes), args.ref_genome, chrom, start, size, args.hal, maf_path)
    system(cmd)
    # phyloFit writes to the file *.mod in whatever the path set by --out-root is
    region_specific_conserved = os.path.join(target.getLocalTempDir(), 'region_specific_conserved')
    cmd = 'phyloFit --init-model {} --scale-only --out-root {} {}'.format(args.conserved_model,
                                                                          region_specific_conserved, maf_path)
    system(cmd)
    region_specific_conserved += '.mod'
    rename_model(region_specific_conserved, args.accelerated_genomes)
    region_specific_accelerated = os.path.join(target.getLocalTempDir(), 'region_specific_accelerated')
    cmd = 'phyloFit --init-model {} --scale-subtree Anc:gain --out-root {} {}'.format(region_specific_conserved,
                                                                                      region_specific_accelerated,
                                                                                      maf_path)
    region_specific_accelerated += '.mod'
    system(cmd)
    region_bed = os.path.join(target.getLocalTempDir(), 'region.bed')
    with open(region_bed, 'w') as outf:
        outf.write('\t'.join([chrom, start, stop]) + '\n')
    cmd = 'phastOdds --output-bed --features {} --background-mods {} --feature-mods {} {} > {}'
    cmd = cmd.format(region_bed, region_specific_conserved, region_specific_accelerated, maf_path, result_path)
    system(cmd)


def cat_results(target, args, paths):
    """
    Concatenates final phastcons output into one gff file.
    """
    fofn = os.path.join(target.getGlobalTempDir(), 'fofn')
    with open(fofn, 'w') as outf:
        for p in paths:
            outf.write(p + '\n')
    tmp_p = os.path.join(target.getGlobalTempDir(), os.path.basename(args.out_bed) + '.tmp')
    cat_cmd = 'cat {} | xargs -n 50 cat > {}'.format(fofn, tmp_p)
    system(cat_cmd)
    os.rename(tmp_p, args.out_bed)


if __name__ == '__main__':
    from phast.run_acceleration_tests import *
    args = parse_args()
    s = Stack(Target.makeTargetFn(extract_maf_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
