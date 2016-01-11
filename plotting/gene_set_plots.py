"""
Produces plots of the protein coding consensus found by consensus.py
"""
import argparse
import os
import cPickle as pickle
import pandas as pd
from collections import OrderedDict, defaultdict
import lib.psl_lib as psl_lib
import lib.sql_lib as sql_lib
import lib.plot_lib as plot_lib
from lib.general_lib import mkdir_p, DefaultOrderedDict, convert_dicts_to_dataframe
import etc.config

__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--refGenome", required=True)
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
        try:
            with open(p) as inf:
                r = pickle.load(inf)
        except IOError:
            continue
        tx_evals[genome] = r["transcript"]
        gene_evals[genome] = r["gene"]
        gene_fail_evals[genome] = r["gene_fail"]
        tx_dup_rate[genome] = r["duplication_rate"]
        tx_counts[genome] = r["tx_counts"]
        gene_counts[genome] = r["gene_counts"]
    return tx_evals, gene_evals, gene_fail_evals, tx_dup_rate, tx_counts, gene_counts


def find_total(data_dict):
    totals = {sum(x.values()) for x in data_dict.itervalues()}
    assert len(totals) == 1, "unequal number of items between genomes"
    return totals.pop()


def transcript_gene_plot(evals, out_path, gencode, mode, biotype):
    results, categories = convert_dicts_to_dataframe(evals, norm=True)
    total = find_total(evals)
    base_title = "Breakdown of {:,} {} {} categorized by consensus finding\nfrom annotation set {}"
    title = base_title.format(total, biotype, mode, gencode)
    out_name = "{}_{}_{}_consensus".format(gencode, biotype, mode)
    palette = etc.config.palette if mode == "genes" or biotype != "protein_coding" else etc.config.triple_palette
    plot_lib.stacked_barplot(results, categories, out_path, out_name, title, color_palette=palette)


def gene_fail_plot(gene_fail_evals, out_path, gencode, biotype):
    results, categories = convert_dicts_to_dataframe(gene_fail_evals)
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
    title = base_title.format(biotype, mode.lower(), gencode)
    out_name = "{}_{}_{}_{}_consensus".format(gencode, biotype, mode, "BinSizes")
    plot_lib.unequal_barplot(results, out_path, out_name, title)


def collapse_evals(tx_evals):
    result = OrderedDict()
    for genome, vals in tx_evals.iteritems():
        tot = sum([y for x, y in vals.iteritems() if x != "NoTransMap"])
        result[genome] = tot
    return result


def biotype_stacked_plot(counter, out_path, gencode, mode):
    results, categories = convert_dicts_to_dataframe(counter)
    if gencode == "GencodePseudoGeneVM7":
        categories = ["\n".join(x.split(" ")) for x in categories]
    base_title = "Biotype breakdown in final {} set derived\nfrom annotation set {}"
    title = base_title.format(mode.lower(), gencode)
    out_name = "{}_{}_{}_consensus".format(gencode, mode, "biotypeStackedPlot")
    plot_lib.stacked_unequal_barplot(results, categories, out_path, out_name, title, 
                                     ylabel="Number of {}s".format(mode.lower()))


def main():
    args = parse_args()
    mkdir_p(args.outDir)
    biotype_tx_counter = DefaultOrderedDict(lambda: defaultdict(int))
    biotype_gene_counter = DefaultOrderedDict(lambda: defaultdict(int))
    gencode_biotype_bin_dict_str = "etc.config.{}".format(args.gencode)
    gencode_biotype_bin_dict = eval(gencode_biotype_bin_dict_str)
    con, cur = sql_lib.attach_databases(args.compAnnPath, mode="reference")
    biotypes = sql_lib.get_all_biotypes(cur, args.refGenome, gene_level=True)
    for biotype in biotypes:
        tx_evals, gene_evals, gene_fail_evals, tx_dup_rate, tx_counts, gene_counts = load_evaluations(args.workDir, args.genomes, biotype)
        if len(tx_evals) == 0:  # a biotype may have nothing
            continue
        if biotype in args.biotypes:
            for (evals, counts), mode in zip(*[[[tx_evals, tx_counts], [gene_evals, gene_counts]], ["transcripts", "genes"]]):
                transcript_gene_plot(evals, args.outDir, args.gencode, mode, biotype)
                size_plot(counts, args.outDir, args.gencode, mode, biotype)
            gene_fail_plot(gene_fail_evals, args.outDir, args.gencode, biotype)
            dup_rate_plot(tx_dup_rate, args.outDir, args.gencode, biotype)
        biotype_bin = gencode_biotype_bin_dict.get(biotype, "Other")
        tx_evals_collapsed = collapse_evals(tx_evals)
        gene_evals_collapsed = collapse_evals(gene_evals)
        for genome, tx_count in tx_evals_collapsed.iteritems():
            biotype_tx_counter[genome][biotype_bin] += tx_count
            biotype_gene_counter[genome][biotype_bin] += gene_evals_collapsed[genome]
    for mode, counter in zip(*[["transcript", "gene"], [biotype_tx_counter, biotype_gene_counter]]):
        biotype_stacked_plot(counter, args.outDir, args.gencode, mode)


if __name__ == "__main__":
    main()
