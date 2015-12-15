import os
import argparse
from collections import Counter, OrderedDict
import numpy as np
import lib.sql_lib as sql_lib
import lib.psl_lib as psl_lib
import lib.plot_lib as plot_lib
from lib.general_lib import mkdir_p
import etc.config


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, nargs="+", required=True, help="genomes in this comparison")
    parser.add_argument("--refGenome", type=str, required=True, help="reference genome")
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--gencode", type=str, required=True, help="current gencode set being analyzed")
    parser.add_argument("--filterChroms", nargs="+", default=["Y", "chrY"], help="chromosomes to ignore")
    return parser.parse_args()


# Hard coded bins used for plots.
paralogy_bins = [0, 1, 2, 3, 4, float('inf')]
coverage_bins = [0, 0.0001, 95.0, 98.0, 99.0, 99.999999, 100.0]
identity_bins = [0, 0.0001, 99.0, 99.5, 99.8, 99.999999, 100.0]


def paralogy(cur, genome):
    """
    Finds the number of paralogous alignments. This is defined as the number of transcript IDs with more than one
    alignment.
    """
    cmd = """SELECT TranscriptId FROM attributes.'{0}'""".format(genome)
    return Counter([x[0] for x in cur.execute(cmd)])


def make_hist(vals, bins, reverse=False, roll=0):
    """
    Makes a histogram out of a value vector given a list of bins. Returns this normalized off the total number.
    Reverse reverses the output relative to bins, roll determines how far to roll the bins around. Useful for putting
    the 0 bin on top.
    """
    raw = np.histogram(vals, bins)[0]
    if reverse is True:
        raw = raw[::-1]
    raw = np.roll(raw, roll)
    norm = raw / (0.01 * len(vals))
    return norm, raw


def get_fail_good_pass_dict(cur, ref_genome, genomes, highest_cov_dict, biotype, filter_chroms):
    """
    Wrapper for sql_lib.get_fail_good_pass_ids that applies to many genomes.
    """
    results = OrderedDict()
    for genome in genomes:
        results[genome] = sql_lib.get_fail_good_pass_ids(cur, ref_genome, genome, biotype, best_cov_only=True,
                                                         filter_chroms=filter_chroms, highest_cov_dict=highest_cov_dict)
    return results


def paralogy_plot(cur, genomes, out_path, biotype, biotype_ids, gencode):
    results = []
    file_name = "{}_{}".format(gencode, "paralogy")
    for g in genomes:
        p = paralogy(cur, g)
        p = [p.get(x, 0) for x in biotype_ids]
        # we roll the list backwards one to put 0 on top
        norm, raw = make_hist(p, paralogy_bins, reverse=False, roll=-1)
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts in {}\nthat have multiple alignments"
    title_string = title_string.format(len(biotype_ids), biotype, gencode)
    legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                    ["= {}".format(paralogy_bins[0])]
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def categorized_plot(cur, highest_cov_dict, genomes, out_path, file_name, biotype, biotype_ids, gencode, query_fn):
    results = []
    for g in genomes:
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        query = query_fn(g, biotype, details=False)
        categorized_ids = sql_lib.get_query_ids(cur, query)
        num_categorized = len({x for x in categorized_ids if x in best_ids and 
                               psl_lib.remove_alignment_number(x) in biotype_ids})
        norm = num_categorized / (0.01 * len(biotype_ids))
        results.append([g, norm, num_categorized])
    title_string = "Proportion of {:,} {} transcripts in biotype {}\ncategorized as {}"
    title_string = title_string.format(len(biotype_ids), biotype, gencode, query_fn.__name__)
    plot_lib.barplot(results, out_path, file_name, title_string, adjust_y=False)


def cat_plot_wrapper(cur, highest_cov_dict, genomes, out_path, biotype, gencode, biotype_ids):
    for query_fn in [etc.config.alignmentErrors, etc.config.assemblyErrors]:
        file_name = "{}_{}".format(gencode, query_fn.__name__)
        categorized_plot(cur, highest_cov_dict, genomes, out_path, file_name, biotype, biotype_ids, gencode,
                         query_fn)


def metrics_plot(highest_cov_dict, bins, genomes, out_path, file_name, biotype, gencode, biotype_ids, analysis):
    results = []
    for g in genomes:
        covs = highest_cov_dict[g]
        vals = [eval(analysis) for tx_id, (aln_id, coverage, identity) in covs.iteritems() if tx_id in biotype_ids]
        vals.extend([0] * (len(biotype_ids) - len(vals)))  # add all of the unmapped transcripts
        norm, raw = make_hist(vals, bins, reverse=True, roll=0)
        results.append([g, norm])
    title_string = "transMap alignment {} breakdown for\n{:,} {} transcripts in biotype {}"
    title_string = title_string.format(analysis, len(biotype_ids), biotype, gencode)
    legend_labels = ["= {0:.1f}%".format(bins[-1])]
    legend_labels.extend(["< {0:.1f}%".format(x) for x in bins[2:-1][::-1]])
    legend_labels.append("= {0:.1f}%".format(bins[0]))
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def cov_ident_wrapper(highest_cov_dict, genomes, out_path,  biotype, gencode, biotype_ids):
    for analysis in ["coverage", "identity"]:
        bins = eval(analysis + "_bins")
        file_name = "{}_{}".format(gencode, analysis)
        metrics_plot(highest_cov_dict, bins, genomes, out_path, file_name, biotype, gencode, biotype_ids, analysis)


def num_good_pass(fail_good_pass_dict, cur, ref_genome, out_path, biotype, gencode, biotype_ids):
    file_name = "{}_num_good_pass".format(gencode)
    results = []
    for genome, (fail_ids, good_specific_ids, pass_ids) in fail_good_pass_dict.iteritems():
        num_no_aln = len(biotype_ids) - sum([len(x) for x in [fail_ids, good_specific_ids, pass_ids]])
        raw = np.array([len(pass_ids), len(good_specific_ids), len(fail_ids), num_no_aln])
        assert all([x >= 0 for x in raw])
        norm = raw / (0.01 * len(biotype_ids))
        results.append([genome, norm])
    title_string = "Proportion of {:,} {} transcripts in biotype {}\ncategorized as Pass/Good/Fail"
    title_string = title_string.format(len(biotype_ids), biotype.replace("_", " "), gencode)
    legend_labels = ["Pass", "Good", "Fail", "NoAln"]
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def num_good_pass_gene_level(fail_good_pass_dict, cur, ref_genome, out_path, biotype, gencode, transcript_gene_map):
    file_name = "{}_num_good_pass_gene_level".format(gencode)
    results = []
    for genome, (fail_ids, good_specific_ids, pass_ids) in fail_good_pass_dict.iteritems():
        pass_genes = {transcript_gene_map[psl_lib.strip_alignment_numbers(x)] for x in pass_ids}
        good_specific_genes = {transcript_gene_map[psl_lib.strip_alignment_numbers(x)] for x in good_specific_ids}
        fail_genes = {transcript_gene_map[psl_lib.strip_alignment_numbers(x)] for x in fail_ids}
        num_genes = len(set(transcript_gene_map.values()))
        num_pass_genes = len(pass_genes)
        num_good_genes = len(good_specific_genes - pass_genes)
        num_fail_genes = len(fail_genes - (good_specific_genes | pass_genes))
        num_no_aln = num_genes - (num_pass_genes + num_good_genes + num_fail_genes)
        raw = np.array([num_pass_genes, num_good_genes, num_fail_genes, num_no_aln])
        assert all([x >= 0 for x in raw])
        norm = raw / (0.01 * num_genes)
        results.append([genome, norm])
    title_string = "Proportion of {:,} {} genes in biotype {}\nwith at least one transcript categorized as Pass/Good/Fail"
    title_string = title_string.format(num_genes, biotype.replace("_", " "), gencode)
    legend_labels = ["Pass", "Good", "Fail", "NoAln"]
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.comparativeAnnotationDir, mode="transMap")
    highest_cov_dict = sql_lib.get_highest_cov_alns(cur, args.genomes, args.filterChroms)
    for biotype in sql_lib.get_all_biotypes(cur, args.refGenome, gene_level=False):
        biotype_ids = sql_lib.get_biotype_ids(cur, args.refGenome, biotype, filter_chroms=args.filterChroms)
        transcript_gene_map = sql_lib.get_transcript_gene_map(cur, args.refGenome, biotype, filter_chroms=args.filterChroms)
        if len(biotype_ids) > 50:  # hardcoded cutoff to avoid issues where this biotype/gencode mix is nearly empty
            fail_good_pass_dict = get_fail_good_pass_dict(cur, args.refGenome, args.genomes, highest_cov_dict, biotype,
                                                          args.filterChroms)
            out_path = os.path.join(args.outDir, "transmap_analysis", biotype)
            mkdir_p(out_path)
            cov_ident_wrapper(highest_cov_dict, args.genomes, out_path, biotype, args.gencode, biotype_ids)
            cat_plot_wrapper(cur, highest_cov_dict, args.genomes, out_path, biotype, args.gencode, biotype_ids)
            paralogy_plot(cur, args.genomes, out_path, biotype, biotype_ids, args.gencode)
            num_good_pass(fail_good_pass_dict, cur, args.refGenome, out_path, biotype, args.gencode, biotype_ids)
            num_good_pass_gene_level(fail_good_pass_dict, cur, args.refGenome, out_path, biotype, args.gencode, 
                                     transcript_gene_map)

if __name__ == "__main__":
    main()