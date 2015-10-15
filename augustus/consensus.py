import argparse
import re
import os
from collections import defaultdict
import lib.sql_lib as sql_lib
import lib.psl_lib as psl_lib
from lib.general_lib import mkdir_p, merge_dicts
import etc.config


__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--binnedTranscriptPath", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--tmGp", required=True)
    parser.add_argument("--filterChroms", nargs="+", default=["Y", "chrY"], help="chromosomes to ignore")
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


def build_data_dict(id_names, id_list, transcript_gene_map):
    data_dict = defaultdict(dict)
    for ids, n in zip(*[id_list, id_names]):
        for aln_id in ids:
            ens_id = psl_lib.strip_alignment_numbers(aln_id)
            gene_id = transcript_gene_map[ens_id]
            if ens_id not in data_dict[gene_id]:
                data_dict[gene_id][ens_id] = defaultdict(list)
            data_dict[gene_id][ens_id][n].append(aln_id)
    return data_dict


def consensus(data_dict, stats, discard_cov_cutoff=0.2, augustus=True):
    for gene, transcript in data_dict.iteritems():
        fail_ids =



def by_coding_biotype(cur, ref_genome, genome, filter_chroms, biotype, transcript_gene_map):
    fail_ids, good_specific_ids, pass_ids = sql_lib.get_fail_good_pass_ids(cur, ref_genome, genome, biotype)
    biotype_ids = sql_lib.get_filtered_biotype_ids(cur, ref_genome, biotype, filter_chroms)
    aug_query = etc.config.augustusEval(genome)
    aug_ids = sql_lib.get_query_ids(cur, aug_query)
    stats = merge_stats(cur, genome)
    id_names = ["fail_ids", "good_specific_ids", "pass_ids", "aug_ids"]
    id_list = [fail_ids, good_specific_ids, pass_ids, aug_ids]
    data_dict = build_data_dict(id_names, id_list, transcript_gene_map)



def by_noncoding_biotype(cur, ref_genome, genome, filter_chroms):
    fail_ids, good_specific_ids, biotype_ids = sql_lib.get_fail_good_biotype_ids(cur, ref_genome, genome, filter_chroms)


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.compAnnPath, mode="augustus")
    transcript_gene_map = sql_lib.get_transcript_gene_map(cur, args.refGenome)
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