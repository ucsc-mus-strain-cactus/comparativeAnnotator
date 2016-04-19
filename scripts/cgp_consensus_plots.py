"""
Produces plots of the protein coding consensus that includes comparative Augustus predictions found by cgp_consensus.py
"""
import os
import cPickle as pickle
from collections import OrderedDict
import pycbio.plotting.plotting as plot_lib
from pycbio.sys.dataOps import munge_nested_dicts_for_plotting, flatten_list_of_lists
from pycbio.sys.fileOps import ensureDir


__author__ = "Ian Fiddes"


def load_evaluations(work_dir, genomes):
    cgp_additions = OrderedDict()
    cgp_replace = OrderedDict()
    new_isoforms = OrderedDict()
    cgp_missing = OrderedDict()
    consensus_stats = OrderedDict()
    removed_cgp = OrderedDict()
    for genome in genomes:
        p = os.path.join(work_dir, genome + ".metrics.pickle")
        with open(p) as inf:
            r = pickle.load(inf)
        cgp_additions[genome] = r["CgpAdditions"]
        cgp_replace[genome] = r["CgpReplace"]
        new_isoforms[genome] = r["NewIsoforms"]
        cgp_missing[genome] = r["CgpAddMissing"]
        consensus_stats[genome] = r["ConsensusStats"]
        removed_cgp[genome] = r['CgpEqualsTMR']
    return cgp_additions, cgp_replace, new_isoforms, cgp_missing, consensus_stats, removed_cgp


def split_side_by_side(data):
    r1, c1 = munge_nested_dicts_for_plotting(OrderedDict([[x, y['CGP']] for x, y in data.iteritems()]))
    r2, c2 = munge_nested_dicts_for_plotting(OrderedDict([[x, y['PacBio']] for x, y in data.iteritems()]))
    results = [r1, r2]
    categories = flatten_list_of_lists([c1, c2])
    return results, categories


def addition_plot(cgp_additions, out_path):
    results, categories = split_side_by_side(cgp_additions)
    title = "Breakdown of the number of new genes/transcripts introduced by Comparative Augustus / PacBio Augustus"
    plot_lib.stacked_side_by_side_unequal_barplot(results, categories, out_path, title, ylabel="Count")


def replace_plot(cgp_replace, out_path):
    results, categories = split_side_by_side(cgp_replace)
    title = "Breakdown of the number of transMap/augustusTMR consensus transcripts replaced by augustusCGP / PacBio Augustus"
    plot_lib.stacked_side_by_side_unequal_barplot(results, categories, out_path, title, ylabel="Count")


def new_isoforms_plot(new_isoforms, out_path):
    results, categories = munge_nested_dicts_for_plotting(new_isoforms, norm=False)
    title = "Breakdown of the number of new isoforms added by Comparative Augustus / PacBio Augustus"
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title)


def missing_plot(cgp_missing, out_path):
    results, categories = split_side_by_side(cgp_missing)
    title = "Breakdown of the number of missing genes/transcripts rescued by Comparative Augustus / PacBio Augustus"
    plot_lib.stacked_side_by_side_unequal_barplot(results, categories, out_path, title)


def consensus_stats_plot(consensus_stats, out_tx, out_gene):
    plot_types = ["Transcript", "Gene"]
    for t, out_path in zip(*[plot_types, [out_tx, out_gene]]):
        data = OrderedDict((x, y[t]) for x, y in consensus_stats.iteritems())
        results, categories = munge_nested_dicts_for_plotting(data, norm=False)
        title = "Breakdown of the origins of the final CGP/TMR/PacBio consensus {} set".format(t.lower())
        ylabel = "Number of {}s".format(t.lower())
        plot_lib.stacked_unequal_barplot(results, categories, out_path, title, ylabel=ylabel)


def match_tmr_plot(removedcgp, out_path):
    results, categories = munge_nested_dicts_for_plotting(removedcgp, norm=False)
    title = 'Number of CGP / PacBio predictions which end up having identical CDS to a TMR and are filtered out'
    plot_lib.side_by_side_unequal_barplot(results, categories, out_path, title)


def generate_consensus_plots(args):
    ensureDir(args.plot_dir)
    cgp_additions, cgp_replace, new_isoforms, cgp_missing, consensus_stats, removed_cgp = load_evaluations(args.metrics_dir,
                                                                                                           args.genomes)
    addition_plot(cgp_additions, args.addition_plot)
    replace_plot(cgp_replace, args.replace_plot)
    new_isoforms_plot(new_isoforms, args.new_isoform_plot)
    missing_plot(cgp_missing, args.missing_plot)
    consensus_stats_plot(consensus_stats, args.consensus_tx_plot, args.consensus_gene_plot)
    match_tmr_plot(removed_cgp, args.consensus_tx_plot)
