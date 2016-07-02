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
from pycbio.sys.mathOps import format_ratio
from comparativeAnnotator.database_queries import get_aln_ids, get_transcript_gene_map
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_db')
    parser.add_argument('ref_genome')
    parser.add_argument('genome')
    parser.add_argument('cgp')
    parser.add_argument('transmap')
    parser.add_argument('--overlap_cutoff', help='Percent of exon overlaps to count', default=0.6, type=float)
    parser.add_argument('--remove-multiple', action='store_true')
    return parser.parse_args()


def create_chrom_dict(tx_dict):
    """
    For all transcripts of a gene on each chromosome, produce a set of exon intervals
    """
    chrom_dict = defaultdict(lambda: defaultdict(set))
    for tx in tx_dict.itervalues():
        chrom_dict[tx.chromosome][tx.name2].update(tx.exon_intervals)
    return chrom_dict


def calculate_overlap(tm_gene_interval, cgp_gene_interval):
    """calculate exon overlaps, returning a total count of matching bases"""
    r = 0
    for tm in tm_gene_interval:
        for cgp in cgp_gene_interval:
            i = tm.intersection(cgp)
            if i is not None:
                r += len(i)
    return r


def find_overlaps(tm_gene_intervals, cgp_gene_interval, min_overlap):
    """
    Compares a cgp record to the transcript intervals for all transcripts on this chromosome.
    Returns only the name of genes which have any transcripts with min_overlap overlap with the cgp_rec.
    """
    r = []
    cgp_size = len(cgp_gene_interval)
    for gene_id, tm_gene_interval in tm_gene_intervals.iteritems():
        overlap = calculate_overlap(tm_gene_interval, cgp_gene_interval)
        percent_overlap = format_ratio(overlap, cgp_size)
        if percent_overlap >= min_overlap:
            r.append([gene_id, percent_overlap])
    return r


def main():
    args = parse_args()
    transmap_dict = get_transcript_dict(args.transmap)
    # pull out best alignment IDs
    best_ids = get_aln_ids(args.ref_genome, args.genome, args.comp_db, best_only=True)
    # filter transMap for large sized transcripts that will mess this up
    transmap_dict = {tx_id: tx for tx_id, tx in transmap_dict.iteritems() if tx_id in best_ids if
                     tx.stop - tx.start <= 3 * 10 ** 6}
    # rename these to the ENSMUSG naming scheme since transMap uses the common names
    transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db)
    for tx_id, tx in transmap_dict.iteritems():
        tx.name2 = transcript_gene_map[strip_alignment_numbers(tx_id)]
    # rename cgp tx's too in case they are not
    cgp_dict = get_transcript_dict(args.cgp)
    cgp_gene_dict = defaultdict(list)
    for tx_id, tx in cgp_dict.iteritems():
        cgp_gene_name = tx.name.split('.')[0]
        tx.name2 = cgp_gene_name
        cgp_gene_dict[cgp_gene_name].append(tx)
    tm_chrom_dict = create_chrom_dict(transmap_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)
    overlap_counts = {}
    for chrom, tm_gene_intervals in tm_chrom_dict.iteritems():
        cgp_gene_intervals = cgp_chrom_dict[chrom]
        for cgp_gene_name, cgp_gene_interval in cgp_gene_intervals.iteritems():
            cgp_recs = cgp_gene_dict[cgp_gene_name]
            r = find_overlaps(tm_gene_intervals, cgp_gene_interval, min_overlap=args.overlap_cutoff)
            overlap_counts[cgp_gene_name] = r  # TODO: make plots out of this
            if args.remove_multiple is True and len(r) > 1:
                continue
            gene_ids, _ = zip(*r)
            for cgp_rec in cgp_recs:
                cgp_rec.name2 = ','.join(gene_ids)
            for cgp_rec in cgp_recs:
                print '\t'.join(map(str, cgp_rec.get_gene_pred()))


if __name__ == '__main__':
    main()
