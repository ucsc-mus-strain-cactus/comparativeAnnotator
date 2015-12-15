"""
Produces plots of the protein coding consensus found by consensus.py
"""
import argparse
import os
import cPickle as pickle
import pandas as pd
from collections import OrderedDict
import lib.psl_lib as psl_lib
import lib.sql_lib as sql_lib
import lib.plot_lib as plot_lib
from lib.general_lib import mkdir_p
import etc.config

__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--gencode", required=True)
    parser.add_argument("--workDir", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--biotypes", default=["protein_coding", "lincRNA", "miRNA", "snRNA", "snoRNA"])
    return parser.parse_args()


def load_evaluations(work_dir, genomes, biotype):
    tx_evals = OrderedDict()
    gene_evals = OrderedDict()
    gene_fail_evals = OrderedDict()
    tx_dup_rate = OrderedDict()
    tx_counts = OrderedDict()
    gene_counts = OrderedDict()
    for genome in genomes:
        p = os.path.join(work_dir, "_".join([genome, biotype]))
        with open(p) as inf:
            r = pickle.load(inf)
        tx_evals[genome] = r["transcript"]
        gene_evals[genome] = r["gene"]
        gene_fail_evals[genome] = r["gene_fail"]
        tx_dup_rate[genome] = r["duplication_rate"]
        tx_counts[genome] = r["tx_counts"]
        gene_counts[genome] = r["gene_counts"]
    return tx_evals, gene_evals, gene_fail_evals, tx_dup_rate, tx_counts, gene_counts


def munge_data(data_dict, norm=False):
    df = pd.DataFrame.from_dict(data_dict)
    categories = data_dict.itervalues().next().keys()
    # hack to reindex rows. from_dict honors only the first layer of ordered dicts
    df = df.reindex(categories)
    if norm is True:
        df = 100 * df.div(df.sum(axis=0), axis=1)
    return [[x[0], x[1].tolist()] for x in df.iteritems()], categories


def find_total(data_dict):
    totals = {sum(x.values()) for x in data_dict.itervalues()}
    assert len(totals) == 1, "unequal number of items between genomes"
    return totals.pop()


def transcript_gene_plot(evals, out_path, gencode, mode, biotype):
    results, categories = munge_data(evals, norm=True)
    total = find_total(evals)
    base_title = "Breakdown of {:,} {} {} categorized by consensus finding\nfrom annotation set {}"
    title = base_title.format(total, biotype, mode, gencode)
    out_name = "{}_{}_{}_consensus".format(gencode, biotype, mode)
    palette = etc.config.palette if mode == "genes" or biotype != "protein_coding" else etc.config.triple_palette
    plot_lib.stacked_barplot(results, categories, out_path, out_name, title, color_palette=palette)


def gene_fail_plot(gene_fail_evals, out_path, gencode, biotype):
    results, categories = munge_data(gene_fail_evals)
    base_title = "Breakdown of {} genes that failed consensus finding\nfrom annotation set {}"
    title = base_title.format(biotype, gencode)
    out_name = "{}_{}_{}_consensus".format(gencode, biotype, "GeneFail")
    plot_lib.stacked_unequal_barplot(results, categories, out_path, out_name, title, ylabel="Number of genes")


def dup_rate_plot(tx_dup_rate, out_path, gencode, biotype):
    results = list(tx_dup_rate.iteritems())
    base_title = "Number of duplicate {} transcripts in consensus before de-duplication\nfrom annotation set {}"
    title = base_title.format(biotype, gencode)
    out_name = "{}_{}_{}_consensus".format(gencode, biotype, "DuplicationRate")
    plot_lib.unequal_barplot(results, out_path, out_name, title)


def size_plot(counts, out_path, gencode, mode, biotype):
    results = list(counts.iteritems())
    base_title = "Number of {} {} in consensus\nfrom annotation set {}"
    title = base_title.format(biotype, mode, gencode)
    out_name = "{}_{}_{}_{}_consensus".format(gencode, biotype, mode, "BinSizes")
    plot_lib.unequal_barplot(results, out_path, out_name, title)


def main():
    args = parse_args()
    mkdir_p(args.outDir)
    for biotype in args.biotypes:
        tx_evals, gene_evals, gene_fail_evals, tx_dup_rate, tx_counts, gene_counts = load_evaluations(args.workDir, args.genomes, biotype)
        for (evals, counts), mode in zip(*[[[tx_evals, tx_counts], [gene_evals, gene_counts]], ["transcripts", "genes"]]):
            transcript_gene_plot(evals, args.outDir, args.gencode, mode, biotype)
            size_plot(counts, args.outDir, args.gencode, mode, biotype)
        gene_fail_plot(gene_fail_evals, args.outDir, args.gencode, biotype)
        dup_rate_plot(tx_dup_rate, args.outDir, args.gencode, biotype)


if __name__ == "__main__":
    main()
