"""
Produces plots of the protein coding consensus that includes comparative Augustus predictions found by cgp_consensus.py
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
from plotting.consensus_plots import munge_data

__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--workDir", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--gencode", required=True)
    return parser.parse_args()


def load_evaluations(work_dir, genomes):
    cgp_additions = OrderedDict()
    cgp_replace = OrderedDict()
    new_isoforms = OrderedDict()
    cgp_missing = OrderedDict()
    cgp_join_genes = OrderedDict()
    for genome in genomes:
        p = os.path.join(work_dir, genome + ".metrics.pickle")
        with open(p) as inf:
            r = pickle.load(inf)
        cgp_additions[genome] = r["CgpAdditions"]
        cgp_replace[genome] = r["CgpReplace"]
        new_isoforms[genome] = r["NewIsoforms"]
        cgp_missing[genome] = r["CgpAddMissing"]
        cgp_join_genes[genome] = r["JoinGeneSupported"]
    return cgp_additions, cgp_replace, new_isoforms, cgp_missing, cgp_join_genes


def addition_plot(cgp_additions, out_path, gencode):
    results, categories = munge_data(cgp_additions, norm=False)
    base_title = ("Breakdown of the number of new genes/transcripts introduced by Comparative Augustus\n"
                  "to the consensus gene set derived from the annotation set {}")
    title = base_title.format(gencode)
    out_name = "{}_{}_cgp_consensus".format(gencode, "gene_addition")
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, out_name, title, ylabel="Count")


def replace_plot(cgp_replace, out_path, gencode):
    results, categories = munge_data(cgp_replace, norm=False)
    base_title = ("Breakdown of the number of transMap/augustusTMR consensus transcripts replaced by augustusCGP\n"
                  "from the consensus gene set derived from the annotation set {}")
    title = base_title.format(gencode)
    out_name = "{}_{}_cgp_consensus".format(gencode, "transcript_replacement")
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, out_name, title, ylabel="Count")


def new_isoforms_plot(new_isoforms, out_path, gencode):
    results = list(new_isoforms.iteritems())
    base_title = ("Breakdown of the number of new isoforms added by Comparative Augustus\n"
                  "to the consensus gene set derived from the annotation set {}")
    title = base_title.format(gencode)
    out_name = "{}_{}_cgp_consensus".format(gencode, "new_isoforms")
    plot_lib.unequal_barplot(results, out_path, out_name, title)


def missing_plot(cgp_missing, out_path, gencode):
    results, categories = munge_data(cgp_missing, norm=False)
    base_title = ("Breakdown of the number of missing genes/transcripts rescued by Comparative Augustus\n"
                  "to the consensus gene set derived from the annotation set {}")
    title = base_title.format(gencode)
    out_name = "{}_{}_cgp_consensus".format(gencode, "missing_genes")
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, out_name, title)


def join_genes_plot(cgp_join_genes, out_path, gencode):
    results, categories = munge_data(cgp_join_genes, norm=False)
    base_title = ("How many CGP consensus transcripts join TMR transcripts in a supported fashion\n"
                  "to the consensus gene set derived from the annotation set {}")
    title = base_title.format(gencode)
    out_name = "{}_{}_cgp_consensus".format(gencode, "join_genes")
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, out_name, title)


def main():
    args = parse_args()
    mkdir_p(args.outDir)
    cgp_additions, cgp_replace, new_isoforms, cgp_missing, cgp_join_genes = load_evaluations(args.workDir, args.genomes)
    addition_plot(cgp_additions, args.outDir, args.gencode)
    replace_plot(cgp_replace, args.outDir, args.gencode)
    new_isoforms_plot(new_isoforms, args.outDir, args.gencode)
    missing_plot(cgp_missing, args.outDir, args.gencode)
    join_genes_plot(cgp_join_genes, args.outDir, args.gencode)


if __name__ == "__main__":
    main()
