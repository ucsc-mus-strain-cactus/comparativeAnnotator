"""
Converts a pileup2fasta output GFF to a PSL representing a fake alignment 

"""

import sys
import os
import argparse
from operator import itemgetter
from itertools import groupby
from collections import defaultdict
from lib.general_lib import tokenize_stream, opener

__author__ = "Ian Fiddes"


class IndelRecord(object):
    """
    Stores indel records from pileup2fasta gff files.
    """
    __slots__ = ["chromosome", "category", "start", "end", "size"]
    def __init__(self, tokens):
        self.chromosome = tokens[0]
        self.category = tokens[2]
        assert self.category in ["DELETION", "INSERTION", "SNP"]
        self.start = int(tokens[3]) - 1
        self.end = int(tokens[4]) - 1
        if self.category == "DELETION":
            self.size = int(tokens[-1].split("=")[-1])
            assert self.size == self.end - self.start
        elif self.category == "INSERTION":
            self.size = len(tokens[-1].split("=")[-1])
        else:
            self.size = 1
    def __repr__(self):
        return "IndelRecord(chromosome={}, start={}, end={}, category={})".format(self.chromosome, self.start, 
                                                                                 self.end, self.category)
    def __len__(self):
        return self.size


class DeletionRecord(IndelRecord):
    def __init__(self, chrom, start, end, size):
        self.chromosome = chrom
        self.start = start
        self.end = end
        self.size = size
        self.category = "DELETION"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refSizes", required=True)
    parser.add_argument("--targetSizes", required=True)
    parser.add_argument("--gff", required=True)
    return parser.parse_args()


def clean_gff_snps(gff_recs, snp_positions):
    """
    These silly GFFs have 1bp deletions listed as SNPs, but the target is a *. It then lines these up sometimes.
    We need to merge these into deletion records.
    """
    intervals = []
    for k, g in groupby(enumerate(snp_positions), lambda (i, (x, y)):i - x):
        intervals.append(map(itemgetter(1), g))
    for interval in intervals:
        positions, chrom = zip(*interval)
        assert len(set(chrom)) == 1
        gff_recs[chrom[0]].append(DeletionRecord(chrom[0], positions[0], positions[0] + len(positions) + 1, len(positions) + 1))


def build_gff_dict(gff_file):
    """
    Builds a sorted list of IndelRecords for each chromosome in the input gff
    Returns a dictionary of chrom: list
    """
    gff_recs = defaultdict(list)
    with open(gff_file) as inf:
        snp_positions = []
        for tokens in tokenize_stream(inf):
            if tokens[2] in ["DELETION", "INSERTION"]:
                rec = IndelRecord(tokens)
                gff_recs[rec.chromosome].append(rec)
            elif tokens[-1][-1] == "*":
                snp_positions.append([int(tokens[3]), tokens[0]])
    # fix the stupid SNP * problem
    clean_gff_snps(gff_recs, snp_positions)
    # just in case for some reason the GFF isn't sorted
    for chrom, val in gff_recs.iteritems():
        gff_recs[chrom] = sorted(val, key=lambda x: x.start)
    return gff_recs


def build_blocks(chrom_recs, q_start, q_end):
    q_num_insert = 0
    q_base_insert = 0
    t_num_insert = 0
    q_starts = [q_start]
    t_starts = [0]
    block_sizes = []
    t_pos = 0
    q_pos = q_start
    for rec in chrom_recs[1:-1]:
        prev_block_size = rec.start - q_starts[-1] + 1
        block_sizes.append(prev_block_size)
        if rec.category == "DELETION":
            t_pos += prev_block_size
            q_pos += rec.end - q_starts[-1]
            q_num_insert += 1
            q_base_insert += rec.size - 1
        elif rec.category == "INSERTION":
            t_pos += prev_block_size + rec.size
            q_pos += prev_block_size
            t_num_insert += 1
        else:
            raise RuntimeError("you shouldn't have gotten here")
        t_starts.append(t_pos)
        q_starts.append(q_pos)
        print rec, len(rec), "q_pos: {}".format(q_pos), "t_pos: {}".format(t_pos), "prev_block: {}".format(prev_block_size)
    block_sizes.append(q_end - q_pos)
    return q_num_insert, q_base_insert, t_num_insert, len(block_sizes), block_sizes, q_starts, t_starts


def gff_to_psl(gff_recs, ref_sizes, target_sizes):
    for chrom, chrom_recs in gff_recs.iteritems():
        t_size = int(target_sizes[chrom])
        q_size = int(ref_sizes[chrom])
        q_start = chrom_recs[0].end
        q_end = chrom_recs[-1].start + 1
        t_start = 0
        t_end = t_size
        vals = build_blocks(chrom_recs, q_start, q_end)
        q_num_insert, q_base_insert, t_num_insert, block_count, block_sizes, q_starts, t_starts = vals
        block_sizes = ",".join(map(str, block_sizes))
        q_starts = ",".join(map(str, q_starts))
        t_starts = ",".join(map(str, t_starts))
        yield [0, 0, 0, 0, q_num_insert, q_base_insert, t_num_insert, 0, "+", chrom, q_size, q_start, q_end,
               chrom, t_size, t_start, t_end, block_count, block_sizes, q_starts, t_starts]


def main():
    args = parse_args()
    ref_sizes = dict(x.split() for x in open(args.refSizes))
    target_sizes = dict(x.split() for x in open(args.targetSizes))
    gff_recs = build_gff_dict(args.gff)
    for x in gff_to_psl(gff_recs, ref_sizes, target_sizes):
        print "\t".join(map(str, x))


if __name__ == "__main__":
    main()
