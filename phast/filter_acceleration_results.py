"""
Filter the LRT results to pick the highest scoring branch, with a minimum score cutoff
"""
import sys
import os
import argparse
from collections import defaultdict

def generate_dict(bed_handle):
    bed_map = defaultdict(list)
    for l in bed_handle:
        chrom, start, stop, branch, score = l.split()
        score = int(score)
        bed_map[(chrom, start, stop)].append([branch, score])
    return bed_map


def find_highest(bed_map, cutoff):
    r = []
    for (chrom, start, stop), recs in bed_map.iteritems():
        recs = sorted(recs, key=lambda (branch, score): -score)
        best_branch, best_score = recs[0]
        if best_score >= cutoff:
            r.append([chrom, start, stop, best_branch, best_score])
    return r


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputbed', type=argparse.FileType('r'))
    parser.add_argument('--cutoff', default=100, type=int)
    args = parser.parse_args()
    bed_map = generate_dict(args.inputbed)
    highest_scoring = find_highest(bed_map, args.cutoff)
    for l in highest_scoring:
        print '\t'.join(map(str, l))
