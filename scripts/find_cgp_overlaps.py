"""
Find CGP overlaps, looking only at transcripts which are tagged as best
"""
import argparse
import sys
import os
import pybedtools
from collections import defaultdict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.intervals import ChromosomeInterval
from comparativeAnnotator.database_queries import get_aln_ids, get_transcript_gene_map
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_db')
    parser.add_argument('ref_genome')
    parser.add_argument('genome')
    parser.add_argument('cgp')
    parser.add_argument('transmap')
    parser.add_argument('--remove-multiple', action='store_true')
    parser.add_argument('--do-not-attempt-to-resolve-multiple', action='store_false')
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


def resolve_multiple(tm_txs, cgp_rec):
    """
    In some cases, we may have multiple assignments. We can try to resolve this by looking at exon overlaps.
    For this, I will use pybedtools.
    """
    cgp = pybedtools.BedTool('\t'.join(map(str, cgp_rec.get_bed())), from_string=True)
    tm = pybedtools.BedTool('\n'.join(['\t'.join(map(str, x.get_bed(name=x.name2))) for x in tm_txs]), from_string=True)
    i = tm.intersect(cgp)
    best = sorted([[x.name, len(x)] for x in i], key=lambda (n, l): -l)[0][0].split('/')[0]
    # because pybedtools sucks
    pybedtools.cleanup()
    return [best]


def evaluate_overlaps(tm_intervals, transmap_dict, cgp_rec, resolve=True, min_overlap=0.2):
    """
    Compares a cgp record to the transcript intervals for all transcripts on this chromosome.
    Returns only the name of genes which have any transcripts with min_overlap overlap with the cgp_rec.
    """
    cgp_interval = ChromosomeInterval(cgp_rec.chromosome, cgp_rec.start, cgp_rec.stop, cgp_rec.strand, cgp_rec.name)
    intersections = [[tm_interval, tm_interval.intersection(cgp_interval)] for tm_interval in tm_intervals]
    overlaps = [[tm_interval, 1.0 * len(intersection) / len(cgp_interval)] for tm_interval, intersection in
                intersections if intersection is not None]
    overlap_names = {tm_interval.name for tm_interval, overlap in overlaps if overlap >= min_overlap}
    if resolve is True and len(overlap_names) > 1:
        tm_txs = [x for x in transmap_dict.itervalues() if x.name2 in overlap_names]
        overlap_names = resolve_multiple(tm_txs, cgp_rec)
    return overlap_names


def main():
    args = parse_args()
    transmap_dict = get_transcript_dict(args.transmap)
    # pull out best alignment IDs
    best_ids = get_aln_ids(args.ref_genome, args.genome, args.comp_db, best_only=True)
    # filter transMap for large sized transcripts that will mess this up
    transmap_dict = {tx_id: tx for tx_id, tx in transmap_dict.iteritems() if tx_id in best_ids if
                     tx.stop - tx.start <= 10 ** 6}
    # rename these to the ENSMUSG naming scheme since transMap uses the common names
    transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db)
    for tx_id, tx in transmap_dict.iteritems():
        tx.name2 = transcript_gene_map[strip_alignment_numbers(tx_id)]
    cgp_dict = get_transcript_dict(args.cgp)
    tm_chrom_dict = create_chrom_dict(transmap_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)
    for chrom, tm_tx_dict in tm_chrom_dict.iteritems():
        tm_intervals = generate_intervals(tm_tx_dict.itervalues())
        for cgp_tx_id, cgp_rec in cgp_chrom_dict[chrom].iteritems():
            r = evaluate_overlaps(tm_intervals, transmap_dict, cgp_rec, resolve=args.do_not_attempt_to_resolve_multiple)
            if args.remove_multiple is True and len(r) > 1:
                continue
            elif len(r) > 0:
                cgp_rec.name2 = ','.join(r)
            else:
                cgp_rec.name2 = cgp_rec.name.split('.')[0]
            print '\t'.join(map(str, cgp_rec.get_gene_pred()))


if __name__ == '__main__':
    main()
