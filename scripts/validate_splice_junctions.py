"""
Validate splice junctions in consensus set with STAR junctions
"""
import argparse
import os
import sys
from bisect import bisect_left
from collections import defaultdict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.bio.intervals import build_intervals_from_bed
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.sys.mathOps import format_ratio
from comparativeAnnotator.comp_lib.name_conversions import aln_id_is_augustus, aln_id_is_transmap


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--consensus_gp', required=True)
    parser.add_argument('--star_junctions', nargs='+', required=True)
    parser.add_argument('--out', required=True, type=argparse.FileType('w'))
    return parser.parse_args()


def binary_search(a, x, lo=0, hi=None):              # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)            # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)                  # find insertion position
    return pos if pos != hi and a[pos] == x else -1  # don't walk off the end


def get_sj_map(sjs):
    r = defaultdict(lambda: defaultdict(list))
    for sj in sjs:
        if sj.strand is True:
            r[sj.chromosome][True].append(sj)
        else:
            r[sj.chromosome][False].append(sj)
    for chrom in r:
        for strand in r[chrom]:
            r[chrom][strand] = sorted(r[chrom][strand], key=lambda x: (x.chromosome, x.start))
    return r


def compare_intervals(tx_intervals, hints_intervals):
    c = 0
    for i in tx_intervals:
        if binary_search(hints_intervals, i) != -1:
            c += 1
    return c


def main():
    args = parse_args()
    tx_dict = get_transcript_dict(args.consensus_gp)
    sjs = set()
    for sj in args.star_junctions:
        sjs.update(build_intervals_from_bed(sj))
    sj_map = get_sj_map(sjs)
    args.out.write('\t'.join(['TxId', 'AlnId', 'Source', 'Ratio', 'NumSupported', 'NumIntrons']) + '\n')
    for tx_id, tx in tx_dict.iteritems():
        r = compare_intervals(tx.intron_intervals, sj_map[tx.chromosome][tx.strand])
        n = len(tx.intron_intervals)
        if aln_id_is_augustus(tx.id):
            s = 'AugustusTMR'
        elif aln_id_is_transmap(tx.id):
            s = 'transMap'
        else:
            s = 'AugustusCGP'
        args.out.write('\t'.join(map(str, [tx_id, tx.id, s, format_ratio(r, n), r, n])) + '\n')


if __name__ == '__main__':
    main()
