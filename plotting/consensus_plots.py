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
    return parser.parse_args()


def load_evaluations(work_dir, genomes):
    tx_evals = OrderedDict()
    gene_evals = OrderedDict()
    gene_fail_evals = OrderedDict()
    tx_dup_rate = OrderedDict()
    for genome in genomes:
        p = os.path.join(work_dir, genome)
        with open(p) as inf:
            r = pickle.load(inf)
        tx_evals[genome] = r["transcript"]
        gene_evals[genome] = r["gene"]
        gene_fail_evals[genome] = r["gene_fail"]
        tx_dup_rate[genome] = r["duplication_rate"]
    return tx_evals, gene_evals, gene_fail_evals, tx_dup_rate


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


def transcript_gene_plot(evals, out_path, gencode, mode):
    results, categories = munge_data(evals, norm=True)
    total = find_total(evals)
    base_title = "Breakdown of {:,} protein-coding {} categorized by consensus finding\nfrom annotation set {}"
    title = base_title.format(total, mode, gencode)
    out_name = "{}_{}_coding_consensus".format(gencode, mode)
    palette = etc.config.palette if mode == "genes" else etc.config.triple_palette
    plot_lib.stacked_barplot(results, categories, out_path, out_name, title, color_palette=palette)


def gene_fail_plot(gene_fail_evals, out_path, gencode):
    results, categories = munge_data(gene_fail_evals)
    base_title = "Breakdown of genes that failed consensus finding\nfrom annotation set {}"
    title = base_title.format(gencode)
    out_name = "{}_{}_coding_consensus".format(gencode, "GeneFail")
    plot_lib.stacked_unequal_barplot(results, categories, out_path, out_name, title, ylabel="Number of genes")


def dup_rate_plot(tx_dup_rate, out_path, gencode):
    results = list(tx_dup_rate.iteritems())
    base_title = "Number of duplicate transcripts in consensus before de-duplication\nfrom annotation set {}"
    title = base_title.format(gencode)
    out_name = "{}_{}_coding_consensus".format(gencode, "DuplicationRate")
    plot_lib.unequal_barplot(results, out_path, out_name, title)


def main():
    args = parse_args()
    mkdir_p(args.outDir)
    tx_evals, gene_evals, gene_fail_evals, tx_dup_rate = load_evaluations(args.workDir, args.genomes)
    for evals, mode in zip(*[[tx_evals, gene_evals], ["transcripts", "genes"]]):
        transcript_gene_plot(evals, args.outDir, args.gencode, mode)
    gene_fail_plot(gene_fail_evals, args.outDir, args.gencode)
    dup_rate_plot(tx_dup_rate, args.outDir, args.gencode)


if __name__ == "__main__":
    main()
