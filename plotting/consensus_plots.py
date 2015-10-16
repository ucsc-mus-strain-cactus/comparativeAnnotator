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


def load_evaluations(work_dir):
    tx_evals = OrderedDict()
    gene_evals = OrderedDict()
    gene_fail_evals = OrderedDict()
    for genome in etc.config.hard_coded_genome_order:
        p = os.path.join(work_dir, genome)
        with open(p) as inf:
            r = pickle.load(inf)
            tx_evals[genome] = r["transcript"]
            gene_evals[genome] = r["gene"]
            gene_fail_evals = r["gene_fail"]
    return tx_evals, gene_evals, gene_fail_evals


def munge_data(data_dict):
    df = pd.DataFrame.from_dict(data_dict)
    categories = data_dict.itervalues().next().keys()
    # hack to reindex rows. from_dict honors only the first layer of ordered dicts
    df = df.reindex(categories)
    return [[x[0], x[1].tolist()] for x in df.iteritems()], categories


def find_total(data_dict):
    totals = {sum(x.values()) for x in data_dict.itervalues()}
    assert len(totals) == 1, "unequal number of items between genomes with equal flag set"
    return totals.pop()


def transcript_gene_plot(evals, out_path, gencode, mode):
    results, categories = munge_data(evals)
    total = find_total(evals)
    base_title = "Breakdown of {:,} Protein-Coding {} Categorized By Consensus Finding\nFrom Annotation Set {}"
    title = base_title.format(total, mode, gencode)
    out_name = "{}_{}_coding_consensus".format(gencode, mode)
    plot_lib.stacked_barplot(results, categories, out_path, out_name, title, color_palette=etc.config.triple_palette)


def gene_fail_plot(gene_fail_evals, out_path, gencode):
    results, categories = munge_data(gene_fail_evals)
    base_title = "Breakdown of genes that failed consensus finding\nFrom Annotation Set {}"
    title = base_title.format(gencode)
    out_name = "{}_{}_coding_consensus".format(gencode, "GeneFail")
    plot_lib.stacked_unequal_barplot(results, categories, out_path, out_name, title)


def main():
    args = parse_args()
    mkdir_p(args.outDir)
    tx_evals, gene_evals, gene_fail_evals = load_evaluations(args.workDir)
    for evals, mode in zip(*[[tx_evals, gene_evals], ["Transcript", "Gene"]]):
        transcript_gene_plot(evals, args.outDir, args.gencode, mode)
    gene_fail_plot(gene_fail_evals, args.outDir, args.gencode)


if __name__ == "__main__":
    main()