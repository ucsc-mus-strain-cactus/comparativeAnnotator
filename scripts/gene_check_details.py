"""
This script builds a extra sqlite3 database that contains the gene-check-details results.
These can be used to compare to the results my own classifiers produced.

Each column is a gene check entry type. If the value is a 1, then gene check found at least
one instance of that problem for this transcript. NULL otherwise.
"""

import os, sys, argparse
from collections import defaultdict

from lib.general_lib import FileType, DirType, FullPaths
import lib.psl_lib
import lib.sqlite_lib as sql_lib    

gene_check_details_ext = ".coding.gene-check-details"

def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str)
    parser.add_argument('--genomes', nargs="+")
    parser.add_argument('--dataDir', type=DirType, action=FullPaths)
    parser.add_argument('--outDb', type=str, default="gene_check_details.db")
    parser.add_argument('--primaryKey', type=str, default="AlignmentID")
    parser.add_argument('--overwriteDb', action="store_true")
    return parser


def parse_details_file(file_path):
    results = defaultdict(set)
    with open(file_path) as f:
        for line in f:
            acc, prob, info, chrom, chromStart, chromEnd = line.rstrip().split("\t")
            results[acc].add(prob)
        return results


def main():
    args = build_parser().parse_args()

    if os.path.exists(args.outDb):
        os.remove(args.outDb)

    for genome in args.genomes:
        results = parse_details_file(os.path.join(args.dataDir, genome + gene_check_details_ext))
        possible_columns = set().union(*results.values())
        columns = [[x, "INTEGER"] for x in possible_columns]
        with sql_lib.ExclusiveSqlConnection(args.outDb) as cur:
            sql_lib.initializeTable(cur, genome, columns, args.primaryKey)
            for aln_id, vals in results.iteritems():
                for val in vals:
                    sql_lib.upsert(cur, genome, args.primaryKey, aln_id, val, "1")


if __name__ == '__main__':
    main()