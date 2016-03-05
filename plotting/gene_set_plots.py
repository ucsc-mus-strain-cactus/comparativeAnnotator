"""
Produces plots of the protein coding consensus found by consensus.py
"""
import os
import cPickle as pickle
from collections import OrderedDict, defaultdict
import comparativeAnnotator.comp_lib.plot_lib as plot_lib
from pycbio.sys.dataOps import convert_dicts_to_dataframe
from pycbio.sys.defaultOrderedDict import DefaultOrderedDict

__author__ = "Ian Fiddes"


def gencode_combinations(geneset):
    """
    hard coded combinations of biotypes to make plots more sane.
    """
    gencode = dict([('snoRNA', 'Small RNAs'), ('snRNA', 'Small RNAs'), ('scaRNA', 'Small RNAs'),
                    ('miRNA', 'Small RNAs'), ('misc_RNA', 'Small RNAs'), ('protein_coding', 'Protein Coding'),
                    ('lincRNA', 'lncRNA'), ('macro_lncRNA', 'lncRNA'), ('bidirectional_promoter_lncrna', 'lncRNA')])
    gencode_pseudo = dict([('processed_pseudogene', 'Processed Pseudogenes'),
                           ('translated_processed_pseudogene', 'Processed Pseudogenes'),
                           ('transcribed_processed_pseudogene', 'Processed Pseudogenes'),
                           ('unprocessed_pseudogene', 'Unprocessed Pseudogenes'),
                           ('translated_unprocessed_pseudogene', 'Unprocessed Pseudogenes'),
                           ('transcribed_unprocessed_pseudogene', 'Unprocessed Pseudogenes')])
    ensembl = dict([('snoRNA', 'Small RNAs'), ('misc_RNA', 'Small RNAs'), ('miRNA', 'miRNA'),
                     ('rRNA', 'rRNA'), ('Mt_rRNA', 'rRNA'), ('protein_coding', 'Protein Coding'),
                     ('Mt_tRNA', 'Small RNAs'), ('snRNA', 'Small RNAs'),
                     ('processed_pseudogene', 'Pseudogene'), ('pseudogene', 'Pseudogene')])
    if geneset.lower() == 'ensembl':
        return ensembl
    elif 'pseudo' in geneset.lower():
        return gencode_pseudo
    else:
        return gencode


def load_evaluations(work_dir, genomes, biotypes):
    for biotype in biotypes:
        d = {'tx_evals': OrderedDict(), 'gene_evals': OrderedDict(), 'gene_fail_evals': OrderedDict(),
             'tx_dup_rate': OrderedDict(), 'tx_counts': OrderedDict(), 'gene_counts': OrderedDict()}
        for genome in genomes:
            p = os.path.join(work_dir, genome + '.pickle')
            try:
                with open(p) as inf:
                    r = pickle.load(inf)
            except IOError:
                continue
            d['tx_evals'][genome] = r[biotype]["transcript"]
            d['gene_evals'][genome] = r[biotype]["gene"]
            d['gene_fail_evals'][genome] = r[biotype]["gene_fail"]
            d['tx_dup_rate'][genome] = r[biotype]["duplication_rate"]
            d['tx_counts'][genome] = r[biotype]["tx_counts"]
            d['gene_counts'][genome] = r[biotype]["gene_counts"]
            yield biotype, d


def find_total(data_dict):
    totals = {sum(x.values()) for x in data_dict.itervalues()}
    assert len(totals) == 1, "unequal number of items between genomes"
    return totals.pop()


def transcript_gene_plot(evals, out_path, mode, biotype):
    results, categories = convert_dicts_to_dataframe(evals, norm=True)
    total = find_total(evals)
    base_title = "Breakdown of {:,} {} {} categorized by consensus finding"
    title = base_title.format(total, biotype, mode)
    palette = plot_lib.palette if mode == "genes" or biotype != "protein_coding" else plot_lib.triple_palette
    plot_lib.stacked_barplot(results, categories, out_path, title, color_palette=palette)


def gene_fail_plot(gene_fail_evals, out_path, biotype):
    results, categories = convert_dicts_to_dataframe(gene_fail_evals)
    base_title = "Breakdown of {} genes that failed consensus finding"
    title = base_title.format(biotype)
    plot_lib.stacked_unequal_barplot(results, categories, out_path, title, ylabel="Number of genes")


def dup_rate_plot(tx_dup_rate, out_path, biotype):
    results = list(tx_dup_rate.iteritems())
    base_title = "Number of duplicate {} transcripts in consensus before de-duplication"
    title = base_title.format(biotype)
    plot_lib.unequal_barplot(results, out_path, title)


def size_plot(counts, out_path, mode, biotype):
    results = list(counts.iteritems())
    base_title = "Number of {} {} in consensus"
    title = base_title.format(biotype, mode.lower())
    plot_lib.unequal_barplot(results, out_path, title)


def collapse_evals(tx_evals):
    result = OrderedDict()
    for genome, vals in tx_evals.iteritems():
        tot = sum([y for x, y in vals.iteritems() if x != "NoTransMap"])
        result[genome] = tot
    return result


def biotype_stacked_plot(counter, out_path, mode):
    results, categories = convert_dicts_to_dataframe(counter)
    base_title = "Biotype breakdown in final {} set"
    title = base_title.format(mode.lower())
    plot_lib.stacked_unequal_barplot(results, categories, out_path, title, ylabel="Number of {}s".format(mode.lower()))


def gene_set_plots(args):
    biotype_map = gencode_combinations(args.gene_set.geneSet)
    biotype_tx_counter = DefaultOrderedDict(lambda: defaultdict(int))
    biotype_gene_counter = DefaultOrderedDict(lambda: defaultdict(int))
    for biotype, d in load_evaluations(args.metrics_dir, args.ordered_target_genomes, args.biotypes):
        plot_cfg = args.biotype_plots[biotype]
        biotype_bin = biotype_map.get(biotype, 'Other')
        transcript_gene_plot(d['tx_evals'], plot_cfg.tx_plot, 'transcripts', biotype_bin)
        transcript_gene_plot(d['gene_evals'], plot_cfg.gene_plot, 'genes', biotype_bin)
        size_plot(d['tx_counts'], plot_cfg.size_plot, 'transcripts', biotype_bin)
        size_plot(d['gene_counts'], plot_cfg.gene_size_plot, 'genes', biotype_bin)
        dup_rate_plot(d['tx_dup_rate'], plot_cfg.dup_rate_plot, biotype_bin)
        gene_fail_plot(d['gene_fail_evals'], plot_cfg.fail_plot, biotype_bin)
        tx_evals_collapsed = collapse_evals(d['tx_evals'])
        gene_evals_collapsed = collapse_evals(d['gene_evals'])
        for genome, tx_count in tx_evals_collapsed.iteritems():
            biotype_tx_counter[genome][biotype_bin] += tx_count
            biotype_gene_counter[genome][biotype_bin] += gene_evals_collapsed[genome]
    for mode, counter, p in [['transcript', biotype_tx_counter, args.transcript_biotype_plot],
                             ['gene', biotype_gene_counter, args.gene_biotype_plot]]:
        biotype_stacked_plot(counter, p, mode)
