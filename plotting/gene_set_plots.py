"""
Produces plots of the protein coding consensus found by consensus.py
"""
import os
import pandas as pd
import cPickle as pickle
from collections import OrderedDict, defaultdict
import pycbio.plotting.plotting as plot_lib
from comparativeAnnotator.generate_gene_set import mode_is_aug
from pycbio.sys.dataOps import munge_nested_dicts_for_plotting, flatten_list_of_lists
from pycbio.sys.defaultOrderedDict import DefaultOrderedDict
__author__ = "Ian Fiddes"


def gencode_combinations(geneset):
    """
    hard coded combinations of biotypes to make plots more sane.
    """
    gencode = dict([('snoRNA', 'small RNAs'), ('snRNA', 'small RNAs'), ('scaRNA', 'small RNAs'),
                    ('miRNA', 'small RNAs'), ('misc_RNA', 'small RNAs'), ('protein_coding', 'protein coding'),
                    ('lincRNA', 'lncRNA'), ('macro_lncRNA', 'lncRNA'), ('bidirectional_promoter_lncrna', 'lncRNA'),
                    ('unitary_pseudogene', 'pseudogene'),
                    ('IG_D_pseudogene', 'pseudogene'),
                    ('IG_C_pseudogene', 'pseudogene'),
                    ('polymorphic_pseudogene', 'pseudogene'),
                    ('TR_J_pseudogene', 'pseudogene'),
                    ('TR_V_pseudogene', 'pseudogene'),
                    ('transcribed_unitary_pseudogene', 'pseudogene'),
                    ('IG_V_pseudogene', 'pseudogene'),
                    ('pseudogene', 'pseudogene'),
                    ('unprocessed_pseudogene', 'pseudogene'),
                    ('transcribed_unprocessed_pseudogene', 'pseudogene'),
                    ('translated_unprocessed_pseudogene', 'pseudogene'),
                    ('transcribed_processed_pseudogene', 'pseudogene'),
                    ('processed_pseudogene', 'pseudogene')])
    gencode_pseudo = dict([('processed_pseudogene', 'Processed pseudogenes'),
                           ('translated_processed_pseudogene', 'Processed pseudogenes'),
                           ('transcribed_processed_pseudogene', 'Processed pseudogenes'),
                           ('unprocessed_pseudogene', 'Unprocessed pseudogenes'),
                           ('translated_unprocessed_pseudogene', 'Unprocessed pseudogenes'),
                           ('transcribed_unprocessed_pseudogene', 'Unprocessed pseudogenes')])
    ensembl = dict([('snoRNA', 'small RNAs'), ('misc_RNA', 'small RNAs'), ('miRNA', 'miRNA'),
                    ('rRNA', 'rRNA'), ('Mt_rRNA', 'rRNA'), ('protein_coding', 'protein coding'),
                    ('Mt_tRNA', 'small RNAs'), ('snRNA', 'small RNAs'),
                    ('processed_pseudogene', 'pseudogene'), ('pseudogene', 'pseudogene')])
    if geneset.lower() == 'ensembl':
        return ensembl
    elif 'pseudo' in geneset.lower():
        return gencode_pseudo
    else:
        return gencode


def load_evaluations(work_dir, genomes, biotypes):
    for biotype in biotypes:
        d = {'tx_evals': OrderedDict(), 'gene_evals': OrderedDict(), 'tx_dup_rate': OrderedDict(),
             'longest_evals': OrderedDict()}
        for genome in genomes:
            p = os.path.join(work_dir, genome + '.pickle')
            try:
                with open(p) as inf:
                    r = pickle.load(inf)
            except IOError:
                continue
            d['tx_evals'][genome] = r[biotype]["transcript"]
            d['gene_evals'][genome] = r[biotype]["gene"]
            d['tx_dup_rate'][genome] = r[biotype]["duplication_rate"]
            d['longest_evals'][genome] =r[biotype]['longest']
        yield biotype, d


def transcript_gene_plot(evals, out_path, mode, biotype, is_consensus):
    results, categories = munge_nested_dicts_for_plotting(evals)
    palette = plot_lib.palette
    if is_consensus:
        base_title = "Breakdown of {} {} categorized by consensus finding"
        if 'protein' in biotype.lower() and mode == 'transcripts':
            palette = plot_lib.paired_palette
    else:
        base_title = "Breakdown of {} {} categorized in transMap gene set"
    title = base_title.format(biotype, mode)
    plot_lib.stacked_unequal_barplot(results, categories, out_path, title, color_palette=palette)


def dup_rate_plot(tx_dup_rate, out_path, biotype, is_consensus):
    results = list(tx_dup_rate.iteritems())
    if is_consensus:
        base_title = "Number of duplicate {} transcripts in consensus geneset before de-duplication"
    else:
        base_title = "Number of duplicate {} transcripts in transMap geneset before de-duplication"
    title = base_title.format(biotype)
    plot_lib.unequal_barplot(results, out_path, title)


def longest_plot(longest_rate, out_path, biotype, is_consensus):
    r1, c1 = munge_nested_dicts_for_plotting(OrderedDict([[x, y['Longest']] for x, y in longest_rate.iteritems()]))
    r2, c2 = munge_nested_dicts_for_plotting(OrderedDict([[x, y['Rescue']] for x, y in longest_rate.iteritems()]))
    results = [r1, r2]
    categories = flatten_list_of_lists([c1, c2])
    if is_consensus:
        base_title = "Rate of longest transcript inclusion/collapsed gene rescue in consensus geneset"
    else:
        base_title = "Rate of longest transcript inclusion/collapsed gene rescue in transMap geneset"
    title = base_title.format(biotype)
    plot_lib.stacked_side_by_side_unequal_barplot(results, categories, out_path, title)


def collapse_evals(tx_evals):
    result = OrderedDict()
    for genome, vals in tx_evals.iteritems():
        tot = sum([y for x, y in vals.iteritems() if x != "NoTransMap"])
        result[genome] = tot
    return result


def biotype_stacked_plot(counter, out_path, mode, is_consensus):
    results, categories = munge_nested_dicts_for_plotting(counter)
    if is_consensus is True:
        base_title = "Biotype breakdown in {} consensus set"
    else:
        base_title = "Biotype breakdown in transMap {} set"
    title = base_title.format(mode.lower())
    plot_lib.stacked_unequal_barplot(results, categories, out_path, title, ylabel="Number of {}s".format(mode.lower()))


def gene_set_plots(args):
    is_consensus = mode_is_aug(args.mode)
    biotype_map = gencode_combinations(args.gene_set.geneSet)
    biotype_tx_counter = DefaultOrderedDict(lambda: defaultdict(int))
    biotype_gene_counter = DefaultOrderedDict(lambda: defaultdict(int))
    for biotype, d in load_evaluations(args.metrics_dir, args.ordered_target_genomes, args.biotypes):
        plot_cfg = args.biotype_plots[biotype]
        biotype_bin = biotype_map.get(biotype, 'other')
        transcript_gene_plot(d['tx_evals'], plot_cfg.tx_plot, 'transcripts', biotype, is_consensus)
        transcript_gene_plot(d['gene_evals'], plot_cfg.gene_plot, 'genes', biotype, is_consensus)
        dup_rate_plot(d['tx_dup_rate'], plot_cfg.dup_rate_plot, biotype, is_consensus)
        longest_plot(d['longest_evals'], plot_cfg.longest_rate_plot, biotype, is_consensus)
        tx_evals_collapsed = collapse_evals(d['tx_evals'])
        gene_evals_collapsed = collapse_evals(d['gene_evals'])
        for genome, tx_count in tx_evals_collapsed.iteritems():
            biotype_tx_counter[genome][biotype_bin] += tx_count
            biotype_gene_counter[genome][biotype_bin] += gene_evals_collapsed[genome]
    for mode, counter, p in [['transcript', biotype_tx_counter, args.transcript_biotype_plot],
                             ['gene', biotype_gene_counter, args.gene_biotype_plot]]:
        biotype_stacked_plot(counter, p, mode, is_consensus)
