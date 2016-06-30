"""
Run phastCons on a series of sub-alignments that represent conserved regions. Run these with the existing
conserved/nonconserved model and also with a modified pair of models which extend the branch lengths for the
nonconserved model.

The final result will be a BED file with the LRT test performed (2 * [alternative likelihood - likelihood]) on each
input site.
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
    parser.add_argument('mafDir', help='Directory containing MAFs split by extract_maf_blocks.')
    parser.add_argument('ref_genome', help='Reference genome.')
    parser.add_argument('conserved_model', help='Original conserved model')
    parser.add_argument('non_conserved_model', help='Original nonconserved model.')
    parser.add_argument('modified_non_conserved_model', help='nonconserved model that was modified.')
    parser.add_argument('out_bed', help='output BED')
    parser.add_argument('--target-coverage', default='0.3',
                        help='target coverage parameter for phastCons. (default: %(default)s)')
    parser.add_argument('--expected-length', default='45',
                        help='expected length parameter for phastCons. (default: %(default)s)')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def phastcons_acceleration_wrapper(target, args):
    """
    Main pipeline wrapper.
    """
    # walk the output dir to find all of the files
    files = []
    for root, dirnames, filenames in os.walk(args.mafDir):
        if filenames:
            files.extend([os.path.join(root, x) for x in filenames if x.startswith('tmp') and
                          os.path.getsize(os.path.join(root, x)) > 0])
    result_dir = target.getGlobalTempDir()
    result_tree = TempFileTree(result_dir)
    assert len(files) > 0
    for f in files:
        genome_loc = os.path.basename(f).split('_', 2)[-1]
        chrom, start, stop = genome_loc.split('_')
        r = result_tree.getTempFile(suffix='_' + genome_loc)
        target.addChildTargetFn(phastcons_accel, args=(args, f, r, chrom, start, stop))
    target.setFollowOnTargetFn(cat_accel, args=(args, result_tree.listFiles()))


def parse_lnl(lnl):
    """parses lnl file"""
    l = open(lnl).next()
    v = l.split(' = ')[-1]
    return float(v)


def phastcons_accel(target, args, f, r, chrom, start, stop):
    """
    Runs phastcons twice, once with the original model pair and once with the new one. Reports a BED of the location
    with the final score, which is (2 * [alternative likelihood - likelihood])
    """
    lnl = os.path.join(target.getLocalTempDir(), 'lnl')
    cmd = 'phastCons {} {},{} --lnl {} --target-coverage {} --expected-length {}'
    cmd = cmd.format(f, args.conserved_model, args.non_conserved_model, lnl, args.target_coverage, args.expected_length)
    system(cmd)
    cons_likelihood = parse_lnl(lnl)
    cmd = 'phastCons {} {},{} --lnl {} --target-coverage {} --expected-length {}'
    cmd = cmd.format(f, args.conserved_model, args.modified_non_conserved_model, lnl, args.target_coverage, args.expected_length)
    system(cmd)
    accel_likelihood = parse_lnl(lnl)
    retval = 2 * (accel_likelihood - cons_likelihood)
    with open(r, 'w') as outf:
        outf.write('\t'.join([chrom, start, stop, str(retval)]) + '\n')


def cat_accel(target, args, paths):
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
    s = Stack(Target.makeTargetFn(phastcons_acceleration_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
