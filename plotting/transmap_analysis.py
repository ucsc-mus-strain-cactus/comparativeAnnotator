import os
import argparse
from collections import Counter
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
coverage_bins = [0, 0.0001, 80.0, 95.0, 99.999999, 100.0]
identity_bins = [0, 0.0001, 99.5, 99.8, 99.999999, 100.0]


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
    Reverse reverse the output relative to bins, roll determines how far to roll the bins around
    """
    raw = np.histogram(vals, bins)[0]
    if reverse is True:
        raw = raw[::-1]
    raw = np.roll(raw, roll)
    norm = raw / (0.01 * len(vals))
    return norm, raw


def paralogy_plot(cur, genome_order, out_path, biotype, biotype_ids, gencode):
    results = []
    file_name = "{}_{}".format(gencode, "paralogy")
    for g in genome_order:
        p = paralogy(cur, g)
        p = [p.get(x, 0) for x in biotype_ids]
        norm, raw = make_hist(p, paralogy_bins, reverse=False, roll=1)
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts in {}\nthat have multiple alignments"
    title_string = title_string.format(len(biotype_ids), biotype, gencode)
    legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                    ["= {}".format(paralogy_bins[0])]
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, biotype_ids, gencode, query_fn):
    results = []
    for g in genome_order:
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        query = query_fn(g, biotype, details=False)
        categorized_ids = sql_lib.get_query_ids(cur, query)
        num_categorized = len({x for x in categorized_ids if x in best_ids})
        norm = num_categorized / (0.01 * len(biotype_ids))
        results.append([g, norm, num_categorized])
    title_string = "Proportion of {:,} {} transcripts in biotype {}\ncategorized as {}"
    title_string = title_string.format(len(biotype_ids), biotype, gencode, query_fn.__name__)
    plot_lib.barplot(results, out_path, file_name, title_string, adjust_y=False)


def cat_plot_wrapper(cur, highest_cov_dict, genome_order, out_path, biotype, gencode, biotype_ids):
    for query_fn in [etc.config.alignmentErrors, etc.config.assemblyErrors]:
        file_name = "{}_{}".format(gencode, query_fn.__name__)
        categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, biotype_ids, gencode,
                         query_fn)


def metrics_plot(highest_cov_dict, bins, genome_order, out_path, file_name, biotype, gencode, biotype_ids, analysis):
    results = []
    for g in genome_order:
        covs = highest_cov_dict[g]
        vals = [eval(analysis) for tx_id, (aln_id, identity, coverage) in covs.iteritems() if tx_id in biotype_ids]
        vals.extend([0] * (len(biotype_ids) - len(vals)))  # add all of the unmapped transcripts
        norm, raw = make_hist(vals, bins, reverse=True, roll=0)
        results.append([g, norm])
    title_string = "transMap alignment {} breakdown for\n{:,} {} transcripts in biotype {}"
    title_string = title_string.format(analysis, len(biotype_ids), biotype, gencode)
    legend_labels = ["= {0:.1f}%".format(bins[-1])]
    legend_labels.extend(["< {0:.1f}%".format(x) for x in bins[2:-1][::-1]])
    legend_labels.append("= {0:.1f}%".format(bins[0]))
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def cov_ident_wrapper(highest_cov_dict, genome_order, out_path,  biotype, gencode, biotype_ids):
    for analysis in ["coverage", "identity"]:
        bins = eval(analysis + "_bins")
        file_name = "{}_{}".format(gencode, analysis)
        metrics_plot(highest_cov_dict, bins, genome_order, out_path, file_name, biotype, gencode, biotype_ids, analysis)


def num_good_pass(highest_cov_dict, cur, genome_order, ref_genome, out_path, biotype, gencode, biotype_ids):
    file_name = "{}_num_good_pass".format(gencode)
    results = []
    for g in genome_order:
        good_query = etc.config.transMapEval(ref_genome, g, biotype, good=True)
        pass_query = etc.config.transMapEval(ref_genome, g, biotype, good=False)
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        good_ids = {x for x in sql_lib.get_query_ids(cur, good_query) if x in best_ids}
        pass_ids = {x for x in sql_lib.get_query_ids(cur, pass_query) if x in best_ids}
        num_no_aln = len(biotype_ids) - len(best_ids)
        num_fail = len(best_ids) - len(good_ids)
        num_good = len(good_ids) - len(pass_ids)
        num_pass = len(pass_ids)
        raw = np.array([num_pass, num_good, num_fail, num_no_aln])
        norm = raw / (0.01 * len(biotype_ids))
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts in biotype {}\ncategorized as Pass/Good/Fail"
    title_string = title_string.format(len(biotype_ids), biotype, gencode)
    legend_labels = ["NoAln", "Fail", "Good", "Pass"]
    plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def get_highest_cov_alns(cur, genomes):
    """
    Dictionary mapping each genome to a dictionary reporting each highest coverage alignment and its metrics
    """
    return {genome: sql_lib.highest_cov_aln(cur, genome) for genome in genomes}


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.comparativeAnnotationDir)
    highest_cov_dict = get_highest_cov_alns(cur, args.genomes)
    # genome_order = plot_lib.find_genome_order(highest_cov_dict, gencode_ids)
    genome_order = etc.config.hard_coded_genome_order
    # we will filter out chromosome Y transcripts for this project
    chr_y_ids = {psl_lib.strip_alignment_numbers(x) for x in sql_lib.get_ids_by_chromosome(cur, args.genome,
                                                                                           args.filterChroms)}
    for biotype in sql_lib.get_all_biotypes(cur, args.refGenome, gene_level=False):
        biotype_ids = sql_lib.filter_biotype_ids(cur, args.refGenome, biotype, chr_y_ids, mode="Transcript")
        if len(biotype_ids) > 50:  # hardcoded cutoff to avoid issues where this biotype/gencode mix is nearly empty
            out_path = os.path.join(args.outDir, biotype)
            mkdir_p(out_path)
            cov_ident_wrapper(highest_cov_dict, genome_order, out_path,biotype, args.gencode, biotype_ids)
            cat_plot_wrapper(cur, highest_cov_dict, genome_order, out_path, biotype, args.gencode, biotype_ids)
            paralogy_plot(cur, genome_order, out_path, biotype, biotype_ids, args.gencode)
            num_good_pass(highest_cov_dict, cur, genome_order, args.refGenome, out_path, biotype,
                          args.gencode, biotype_ids)


if __name__ == "__main__":
    main()