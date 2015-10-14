import argparse
import re
import os
import itertools
from collections import defaultdict
from plotting.plot_functions import get_all_biotypes, get_gene_map, get_gene_biotype_map, get_transcript_biotype_map, \
                                    gp_chrom_filter, get_all_ids, get_reverse_name_map
import cPickle as pickle
import lib.sql_lib as sql_lib
import lib.psl_lib as psl_lib
from lib.general_lib import mkdir_p, merge_dicts


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--binnedTranscriptPath", required=True)
    parser.add_argument("--attributePath", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--tmGp", required=True)
    parser.add_argument("--compGp", required=True)
    parser.add_argument("--basicGp", required=True)
    return parser.parse_args()


def merge_stats(cur, genome):
    """
    Adapter function to return the combination of the stats dicts produced for augustus and transMap
    """
    tm_stats = sql_lib.get_stats(cur, genome, mode="transMap")
    aug_stats = sql_lib.get_stats(cur, genome, mode="augustus")
    return merge_dicts([tm_stats, aug_stats])


def find_best_aln(stats, sig_fig=4):
    """
    Takes a list of stats for transcripts and finds all transcripts which have the highest percent identity
    """
    s = sorted(stats, key=lambda x: -x[1])
    best_ident = round(s[0][1], sig_fig)
    return {name for name, aln_id, aln_cov in s if round(aln_id, sig_fig) == best_ident}


def get_good_pass(cur, genome, biotype):
    good_query = etc.config.transMapEval(genome, biotype, good=True)
    pass_query = etc.config.transMapEval(genome, biotype, good=False)
    best_ids = set(zip(*highest_cov_dict[genome].itervalues())[0])
    good_ids = {x for x in sql_lib.get_query_ids(cur, good_query) if strip_alignment_numbers(x) in filter_set
                and x in best_ids}
    pass_ids = {x for x in sql_lib.get_query_ids(cur, pass_query) if strip_alignment_numbers(x) in filter_set
                    and x in best_ids}

def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.compAnnPath, mode="augustus")
    biotypes = get_all_biotypes(args.attributePath)
    gene_map = get_gene_map(args.attributePath)  # maps each transcript to its parent gene
    gene_biotype_map = get_gene_biotype_map(args.attributePath)  # maps each gene to its biotype
    transcript_biotype_map = get_transcript_biotype_map(args.attributePath)  # maps each transcript to its biotype
    chr_y_ids = gp_chrom_filter(args.compGp)  # we want to filter out chrY
    consensus_base_path = os.path.join(args.outDir, "geneSets")
    mkdir_p(consensus_base_path)
    gps = load_gps([args.tmGp, args.augGp])  # load all Augustus and transMap transcripts into one big dict


if __name__ == "__main__":
    main()