"""
Produce two sets of biotype names - one list are those we could not map at all, one list are those that failed
gene set production.
"""
import os
import sys
import argparse
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.database_queries import get_aln_ids, get_ref_ids
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('ref_genome')
    parser.add_argument('db_path')
    parser.add_argument('consensus_gp')
    parser.add_argument('out_list', type=argparse.FileType('w'))
    parser.add_argument('--biotype', default='lincRNA')
    parser.add_argument('--filterChroms', default=['Y', 'chrY'])
    return parser.parse_args()


def main():
    args = parse_args()
    ref_ids = get_ref_ids(args.ref_genome, args.db_path, args.biotype, args.filterChroms)
    aln_ids = get_aln_ids(args.ref_genome, args.genome, args.db_path, args.biotype, best_cov_only=True)
    aln_ids = {strip_alignment_numbers(x) for x in aln_ids}
    gene_set_ids = set(get_transcript_dict(args.consensus_gp).keys())
    assert len(ref_ids & gene_set_ids) == len(gene_set_ids)
    no_tm = ref_ids - aln_ids
    fail_ids = aln_ids - gene_set_ids
    for aln_id in no_tm:
        l = [aln_id, 'NoTransMap']
        args.out_list.write('\t'.join(l) + '\n')
    for aln_id in fail_ids:
        l = [aln_id, 'FailConsensus']
        args.out_list.write('\t'.join(l) + '\n')


if __name__ == '__main__':
    main()
