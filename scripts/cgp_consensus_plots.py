"""
Produces plots of the protein coding consensus that includes comparative Augustus predictions found by cgp_consensus.py
"""
import os
import cPickle as pickle
from collections import OrderedDict
import pycbio.plotting.plotting as plot_lib
from pycbio.sys.dataOps import munge_nested_dicts_for_plotting
from pycbio.sys.fileOps import ensureDir


__author__ = "Ian Fiddes"


def load_evaluations(work_dir, genomes):
    cgp_additions = OrderedDict()
    cgp_replace = OrderedDict()
    new_isoforms = OrderedDict()
    cgp_missing = OrderedDict()
    cgp_join_genes = OrderedDict()
    consensus_stats = OrderedDict()
    for genome in genomes:
        p = os.path.join(work_dir, genome + ".metrics.pickle")
        with open(p) as inf:
            r = pickle.load(inf)
        cgp_additions[genome] = r["CgpAdditions"]
        cgp_replace[genome] = r["CgpReplace"]
        new_isoforms[genome] = r["NewIsoforms"]
        cgp_missing[genome] = r["CgpAddMissing"]
        cgp_join_genes[genome] = r["JoinGeneSupported"]
        consensus_stats[genome] = r["ConsensusStats"]
    return cgp_additions, cgp_replace, new_isoforms, cgp_missing, cgp_join_genes, consensus_stats


def addition_plot(cgp_additions, out_path):
    results, categories = munge_nested_dicts_for_plotting(cgp_additions, norm=False)
    title = "Breakdown of the number of new genes/transcripts introduced by Comparative Augustus"
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title, ylabel="Count")


def replace_plot(cgp_replace, out_path):
    results, categories = munge_nested_dicts_for_plotting(cgp_replace, norm=False)
    title = "Breakdown of the number of transMap/augustusTMR consensus transcripts replaced by augustusCGP"
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title, ylabel="Count")


def new_isoforms_plot(new_isoforms, out_path):
    results = list(new_isoforms.iteritems())
    title = "Breakdown of the number of new isoforms added by Comparative Augustus"
    plot_lib.unequal_barplot(results, out_path, title)


def missing_plot(cgp_missing, out_path,):
    results, categories = munge_nested_dicts_for_plotting(cgp_missing, norm=False)
    title = "Breakdown of the number of missing genes/transcripts rescued by Comparative Augustus"
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title)


def join_genes_plot(cgp_join_genes, out_path):
    results, categories = munge_nested_dicts_for_plotting(cgp_join_genes, norm=False)
    title = "How many CGP consensus transcripts join TMR transcripts in a supported fashion"
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title)


def consensus_stats_plot(consensus_stats, out_path):
    # make this pandas-multi-indexable. TODO: do everything in pandas.
    #reform = {(outer_key, inner_key): values for outer_key, inner_dict in consensus_stats.iteritems() for 
    #          inner_key, values in inner_dict.iteritems()}
    categories = ["Transcript", "Gene"]
    for cat in categories:
        data = OrderedDict((x, y[cat]) for x, y in consensus_stats.iteritems())
        results, categories = munge_nested_dicts_for_plotting(data, norm=False)
        base_title = "Breakdown of the origins of the final CGP/TMR consensus set"
        plot_lib.stacked_unequal_barplot(results, categories, out_path, ylabel="Number of {}".format(cat))


def generate_consensus_plots(args):
    ensureDir(args.plot_dir)
    cgp_additions, cgp_replace, new_isoforms, cgp_missing, cgp_join_genes, consensus_stats = load_evaluations(args.args_holder,
                                                                                                              args.genomes)
    addition_plot(cgp_additions, args.addition_plot)
    replace_plot(cgp_replace, args.replace_plot)
    new_isoforms_plot(new_isoforms, args.new_isoform_plot)
    missing_plot(cgp_missing, args.missing_plot)
    join_genes_plot(cgp_join_genes, args.join_genes_plot)
    consensus_stats_plot(consensus_stats, args.consensus_stats_plot)
