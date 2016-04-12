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
from comparativeAnnotator.database_queries import get_gene_biotype_map


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_db')
    parser.add_argument('ref_genome')
    parser.add_argument('cgp')
    parser.add_argument('tmr')
    return parser.parse_args()


def create_chrom_dict(tx_dict):
    chrom_dict = defaultdict(dict)
    for tx_id, tx in tx_dict.iteritems():
        chrom_dict[tx.chromosome][tx_id] = tx
    return chrom_dict


def generate_intervals(tx_iter):
    return [ChromosomeInterval(x.chromosome, x.start, x.stop, x.strand, x.name2) for x in tx_iter]


def evaluate_overlaps(tm_intervals, cgp_rec):
    cgp_interval = ChromosomeInterval(cgp_rec.chromosome, cgp_rec.start, cgp_rec.stop, cgp_rec.strand, cgp_rec.name)
    return {x.name for x in tm_intervals if x.overlap(cgp_interval, stranded=True)}


def main():
    args = parse_args()
    gene_biotype_map = get_gene_biotype_map(args.ref_genome, args.comp_db)
    not_ok_txs = {gene_id for gene_id, biotype in gene_biotype_map.iteritems() if biotype not in ['protein_coding',
                                                                                                  'processed_transcript']}
    tm_dict = get_transcript_dict(args.tmr)
    cgp_dict = get_transcript_dict(args.cgp)
    tm_chrom_dict = create_chrom_dict(tm_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)
    interesting = defaultdict(list)
    for chrom, tm_tx_dict in tm_chrom_dict.iteritems():
        tm_intervals = generate_intervals(tm_tx_dict.itervalues())
        for cgp_tx_id, cgp_rec in cgp_chrom_dict[chrom].iteritems():
            r = evaluate_overlaps(tm_intervals, cgp_rec)
            r -= not_ok_txs
            interesting[len(r)].append(cgp_rec)
            if len(r) > 0:
                cgp_rec.name2 = ','.join(r)
            else:
                cgp_rec.name2 = cgp_rec.name.split('.')[0]
            print '\t'.join(map(str, cgp_rec.get_gene_pred()))


if __name__ == '__main__':
    main()
