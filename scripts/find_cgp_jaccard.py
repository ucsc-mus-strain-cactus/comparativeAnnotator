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
from pycbio.sys.fileOps import TemporaryFilePath
from pycbio.sys.procOps import callProcLines
from comparativeAnnotator.database_queries import get_aln_ids, get_transcript_gene_map, get_gene_biotype_map
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_db')
    parser.add_argument('ref_genome')
    parser.add_argument('genome')
    parser.add_argument('cgp')
    parser.add_argument('transmap')
    res = parser.add_mutually_exclusive_group()
    res.add_argument('--remove-multiple', action='store_true')
    res.add_argument('--resolve-multiple', action='store_true')
    return parser.parse_args()


def create_chrom_dict(tx_dict):
    """
    Split up a dictionary of Transcript objects by chromosome, converting each object into a BED record
    """
    chrom_dict = defaultdict(dict)
    for tx_id, tx in tx_dict.iteritems():
        chrom_dict[tx.chromosome][tx_id] = [tx, '\t'.join(map(str, tx.get_bed())) + '\n']
    return chrom_dict


def filter_tm_txs(cgp_tx, tm_tx_dict):
    """reduce transcripts to those who intersect the cgp_tx"""
    cgp_interval = ChromosomeInterval(cgp_tx.chromosome, cgp_tx.start, cgp_tx.stop, cgp_tx.strand)
    r = []
    for tx, bed_rec in tm_tx_dict.itervalues():
        tx_interval = ChromosomeInterval(tx.chromosome, tx.start, tx.stop, tx.strand)
        if tx_interval.intersection(cgp_interval) is not None:
            r.append([tx, bed_rec])
    return r


def resolve_multiple_genes(results, gene_biotype_map, min_distance=0.2):
    """
    Resolve multiple assignments based on the following rules:
    If a CGP has multiple assignments, see if only 1 is not a pseudogene. If so, assign it.
    If the difference in Jaccard scores is >=min_distance, assign it to the higher score.
    If neither of these are satisfiable, discard the transcript.
    """
    filtered = {gene_id: score for gene_id, score in results.iteritems() if gene_biotype_map[gene_id] == 'protein_coding'}
    if len(filtered) == 0:
        return None
    elif len(filtered) == 1:
        return filtered.keys()[0]
    scores = results.values()
    high_score = max(scores)
    if all(high_score - x >= min_distance for x in scores if x != high_score):
        best = sorted(results.iteritems(), key=lambda (gene_id, score): score)[-1]
        results = best[0]
    else:
        results = None
    return results


def calculate_jaccard(cgp_bed_rec, filtered_tm_txs, gene_biotype_map):
    """calculates jaccard distance. the pybedtools wrapper can't do stranded"""
    results = defaultdict(float)
    with TemporaryFilePath() as cgp, TemporaryFilePath() as tm:
        for tm_tx, tm_bed_rec in filtered_tm_txs:
            with open(cgp, 'w') as outf:
                outf.write(cgp_bed_rec)
            with open(tm, 'w') as outf:
                outf.write(tm_bed_rec)
            cmd = ['bedtools', 'jaccard', '-s', '-a', cgp, '-b', tm]
            r = callProcLines(cmd)
            j = float(r[-1].split()[-2])
            results[tm_tx.name2] = max(results[tm_tx.name2], j)
    results = resolve_multiple_genes(results, gene_biotype_map)
    return results


def find_overlapping_genes(filtered_tm_txs):
    """
    Determine if transMap overlaps lead to only one gene assocation
    """
    gene_ids = {tx.name2 for tx, bed_rec in filtered_tm_txs}
    return gene_ids


def main():
    args = parse_args()
    transmap_dict = get_transcript_dict(args.transmap)
    # pull out best alignment IDs
    best_ids = get_aln_ids(args.ref_genome, args.genome, args.comp_db, best_only=True)
    transmap_dict = {tx_id: tx for tx_id, tx in transmap_dict.iteritems() if tx_id in best_ids}
    # rename these to the ENSMUSG naming scheme since transMap uses the common names
    transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db)
    gene_biotype_map = get_gene_biotype_map(args.ref_genome, args.comp_db)
    for tx_id, tx in transmap_dict.iteritems():
        tx.name2 = transcript_gene_map[strip_alignment_numbers(tx_id)]
    cgp_dict = get_transcript_dict(args.cgp)
    tm_chrom_dict = create_chrom_dict(transmap_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)
    for chrom, tm_tx_dict in tm_chrom_dict.iteritems():
        for cgp_tx_id, (cgp_tx, cgp_bed_rec) in cgp_chrom_dict[chrom].iteritems():
            filtered_tm_txs = filter_tm_txs(cgp_tx, tm_tx_dict)
            gene_ids = find_overlapping_genes(filtered_tm_txs)
            if len(gene_ids) == 0:
                cgp_tx.name2 = cgp_tx.name.split('.')[0]
            elif len(gene_ids) == 1:
                cgp_tx.name2 = list(gene_ids)[0]
            elif args.remove_multiple is True:
                continue
            elif args.resolve_multiple is True:
                gene_id = calculate_jaccard(cgp_bed_rec, filtered_tm_txs, gene_biotype_map)
                if gene_id is None:  # we cannot resolve this transcript
                    continue
                cgp_tx.name2 = gene_id
            else:
                cgp_tx.name2 = ','.join(gene_ids)
            print '\t'.join(map(str, cgp_tx.get_gene_pred()))


if __name__ == '__main__':
    main()
