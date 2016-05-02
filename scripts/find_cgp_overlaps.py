"""
Find CGP overlaps, looking only at transcripts which are tagged as best
"""
import argparse
import sys
import os
from collections import defaultdict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.intervals import ChromosomeInterval
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers
from comparativeAnnotator.database_queries import get_transcript_biotype_map


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_db')
    parser.add_argument('ref_genome')
    parser.add_argument('cgp')
    parser.add_argument('tmr')
    parser.add_argument('--remove-multiple', action='store_true')
    return parser.parse_args()


def create_chrom_dict(tx_dict):
    """
    Creates a dictionary mapping chromosome names to all transcripts on it
    """
    chrom_dict = defaultdict(dict)
    for tx_id, tx in tx_dict.iteritems():
        chrom_dict[tx.chromosome][tx_id] = tx
    return chrom_dict


def generate_intervals(tx_iter):
    """
    Uses the gene information for a iterable of transcript objects to produce a set of intervals.
    """
    return {ChromosomeInterval(tx.chromosome, tx.start, tx.stop, tx.strand, tx.name2) for tx in tx_iter}


def evaluate_overlaps(tm_intervals, cgp_rec, min_overlap=0.75):
    """
    Compares a cgp record to the transcript intervals for all transcripts on this chromosome.
    Returns only the name of genes which have any transcripts with min_overlap overlap with the cgp_rec.
    """
    cgp_interval = ChromosomeInterval(cgp_rec.chromosome, cgp_rec.start, cgp_rec.stop, cgp_rec.strand, cgp_rec.name)
    intersections = [[tm_interval, tm_interval.intersection(cgp_interval)] for tm_interval in tm_intervals]
    overlaps = [[tm_interval, 1.0 * len(intersection) / len(tm_interval)] for tm_interval, intersection in
                intersections if intersection is not None]
    return {tm_interval.name for tm_interval, overlap in overlaps if overlap >= min_overlap}


def main():
    args = parse_args()
    transcript_biotype_map = get_transcript_biotype_map(args.ref_genome, args.comp_db)
    tm_dict = get_transcript_dict(args.tmr)
    tm_dict = {tx_id: tx for tx_id, tx in tm_dict.iteritems() if transcript_biotype_map[strip_alignment_numbers(tx_id)] == 'protein_coding'}
    cgp_dict = get_transcript_dict(args.cgp)
    tm_chrom_dict = create_chrom_dict(tm_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)
    for chrom, tm_tx_dict in tm_chrom_dict.iteritems():
        tm_intervals = generate_intervals(tm_tx_dict.itervalues())
        for cgp_tx_id, cgp_rec in cgp_chrom_dict[chrom].iteritems():
            r = evaluate_overlaps(tm_intervals, cgp_rec)
            if args.remove_multiple is True and len(r) > 1:
                continue
            elif len(r) > 0:
                cgp_rec.name2 = ','.join(r)
            else:
                cgp_rec.name2 = cgp_rec.name.split('.')[0]
            print '\t'.join(map(str, cgp_rec.get_gene_pred()))


if __name__ == '__main__':
    main()
