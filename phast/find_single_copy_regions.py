"""
Run acceleration tests on each region in the conserved_bed.
"""

import os
import argparse
import itertools
from pyfasta import Fasta
from ete3 import Tree
from collections import OrderedDict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, TempFileTree, popenCatch


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hal', help='HAL alignment file.', required=True)
    parser.add_argument('--ref_genome', help='Reference genome.', required=True)
    parser.add_argument('--conserved_bed', help='Conserved BED file.', required=True)
    parser.add_argument('--out_bed', help='output BED', required=True)
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def format_ratio(numerator, denominator):
    """
    Convenience function that converts two numbers, integer or no, to a ratio
    """
    if denominator == 0:
        return float("nan")
    return float(numerator) / denominator


def grouper(iterable, size):
    """
    http://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
    """
    it = iter(iterable)
    chunk = tuple(itertools.islice(it, size))
    while chunk:
        yield chunk
        chunk = tuple(itertools.islice(it, size))


def single_copy_wrapper(target, args):
    """
    Main pipeline wrapper. Runs halSingleCopyRegionsExtract once for each region in the conserved_bed file.
    """
    bed_recs = [x.split()[:3] for x in open(args.conserved_bed)]
    result_dir = target.getGlobalTempDir()
    result_tree = TempFileTree(result_dir)
    for chunk in grouper(bed_recs, 10):
        result_path = result_tree.getTempFile()
        target.addChildTargetFn(find_single_copy, args=(args, chunk, result_path))
    target.setFollowOnTargetFn(cat_results, args=(args, result_tree.listFiles()))


def find_single_copy(target, args, chunk, result_path):
    """
    Score each region for percent single copyness
    """
    with open(result_path, 'w') as outf:
        for chrom, start, stop in chunk:
            start = int(start)
            stop = int(stop)
            length = stop - start
            cmd = 'halSingleCopyRegionsExtract {} {} --refSequence {} --start {} --length {}'
            cmd = cmd.format(args.hal, args.ref_genome, chrom, start, length)
            r = popenCatch(cmd)
            r = r.split('\n')[:-1]
            tot = 0
            for l in r:
                l = l.split()
                tot += int(l[-1]) - int(l[-2])
            outf.write('\t'.join(map(str, [chrom, start, stop, format_ratio(tot, length)])) + '\n')


def cat_results(target, args, paths):
    """
    Concatenates final scores output into one bed file.
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
    from phast.find_single_copy_regions import *
    args = parse_args()
    s = Stack(Target.makeTargetFn(single_copy_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
