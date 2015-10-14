import argparse
from collections import Counter
import numpy as np
import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
from lib.general_lib import mkdir_p
from plotting.plot_functions import *
import etc.config


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, nargs="+", required=True, help="genomes in this comparison")
    parser.add_argument("--refGenome", type=str, required=True, help="reference genome")
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--gencode", type=str, required=True, help="current gencode set being analyzed")
    return parser.parse_args()


# Hard coded bins used for plots.
paralogy_bins = [0, 1, 2, 3, 4, float('inf')]
identity_bins = [0, 0.0001, 0.995, 0.998, 0.99999999, 1.0]
coverage_bins = [0, 0.0001, 0.8, 0.95, 0.99999999, 1.0]


def paralogy(cur, genome):
    """
    Finds the number of paralogous alignments. This is defined as the number of transcript IDs with more than one
    alignment.
    """
    cmd = """SELECT TranscriptId FROM attributes.'{0}'""".format(genome)
    return Counter([x[0] for x in cur.execute(cmd).fetchall()])


def find_genome_order(highest_cov_dict, filter_set):
    """
    Finds the genome order that will be used by all plots. This is deprecated in favor of using a hard coded order
    provided by Joel.
    """
    num_cov = {}
    for g, covs in highest_cov_dict.iteritems():
        num_cov[g] = 1.0 * len({tx_id for tx_id, (aln_id, cov, ident) in covs.iteritems() if tx_id in filter_set})
        num_cov[g] /= len(filter_set)
    order = sorted(num_cov.iteritems(), key=lambda x: -x[1])
    return zip(*order)[0]


def make_hist(vals, bins, total, reverse=False, roll=True):
    raw = np.histogram(vals, bins)[0]
    if reverse is True:
        raw = raw[::-1]
    if roll is True:
        raw = np.roll(raw, 1)
    norm = raw / (0.01 * total)
    return norm, raw


def paralogy_plot(cur, genome_order, ref_genome, out_path, base_file_name, biotype, gencode):
    results = []
    file_name = "{}_{}".format(base_file_name, "paralogy")
    biotype_ids = sql_lib.get_biotype_ids(cur, ref_genome, biotype)
    for g in genome_order:
        p = paralogy(cur, g)
        p = [p.get(x, 0) for x in biotype_ids]
        norm, raw = make_hist(p, paralogy_bins, len(biotype_ids), reverse=False, roll=True)
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts in {}\nthat have multiple alignments".format(len(biotype_ids), 
                                                                                                   biotype, gencode)
    legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                     ["= {}".format(paralogy_bins[0])] 
    stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, gencode, filter_set, query_fn):
    results = []
    for g in genome_order:
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        r = number_categorized(cur, g, query_fn, biotype=biotype)
        raw = len({x for x in r if strip_alignment_numbers(x) in filter_set and x in best_ids})
        norm = raw / (0.01 * len(filter_set))
        results.append([g, norm, raw])
    title_string = "Proportion of {:,} {} transcripts in {}\ncategorized as {}".format(len(filter_set), biotype, 
                                                                                       gencode, query_fn.__name__)
    barplot(results, out_path, file_name, title_string, adjust_y=False)


def cat_plot_wrapper(cur, highest_cov_dict, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    for query_fn in [etc.config.alignmentErrors, etc.config.assemblyErrors]:
        file_name = "{}_{}".format(base_file_name, query_fn.__name__)
        categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, gencode, filter_set,
                         query_fn)


def metrics_plot(highest_cov_dict, bins, genome_order, out_path, file_name, biotype, gencode, filter_set, analysis):
    results = []
    for g in genome_order:
        mets = highest_cov_dict[g]
        mets = [eval(analysis) for tx_id, (aln_id, identity, coverage) in mets.iteritems() if tx_id in filter_set]
        mets.extend([0] * (len(filter_set) - len(mets)))
        g, norm, raw = make_hist(mets, bins, len(filter_set), g, reverse=True)
        results.append([g, norm])
    title_string = "transMap alignment {} breakdown for\n{:,} {} transcripts in {}".format(analysis, len(filter_set),
                                                                                                 biotype, gencode)
    legend_labels = ["= {0:.1f}%".format(100 * bins[-1])]
    legend_labels.extend(["< {0:.1f}%".format(100 * x) for x in bins[2:-1][::-1]])
    legend_labels.append("= {0:.1f}%".format(100 * bins[0]))
    stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def cov_ident_wrapper(highest_cov_dict, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    for analysis in ["coverage", "identity"]:
        bins = eval(analysis + "_bins")
        file_name = "{}_{}".format(base_file_name, analysis)
        metrics_plot(highest_cov_dict, bins, genome_order, out_path, file_name, biotype, gencode, filter_set, analysis)


def num_good_pass(highest_cov_dict, cur, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    file_name = "{}_numGoodPass".format(base_file_name)
    coding = True if biotype == "protein_coding" else False
    results = []
    for genome in genome_order:
        good_query = etc.config.transMapEval(genome, coding, good=True)
        pass_query = etc.config.transMapEval(genome, coding, good=False)
        best_ids = set(zip(*highest_cov_dict[genome].itervalues())[0])
        good_ids = {x for x in sql_lib.get_query_ids(cur, good_query) if strip_alignment_numbers(x) in filter_set
                    and x in best_ids}
        pass_ids = {x for x in sql_lib.get_query_ids(cur, pass_query) if strip_alignment_numbers(x) in filter_set
                    and x in best_ids}
        num_fail = len(filter_set) - len(best_ids)
        num_pass = len(pass_ids - good_ids)
        num_good = len(good_ids)
        raw = len({x for x in tm_ok if strip_alignment_numbers(x) in filter_set and x in best_ids})
        norm = raw / (0.01 * len(filter_set))
        results.append([genome, norm, raw])
    title_string = "Proportion of {:,} {} transcripts in {}\ncategorized as OK".format(len(filter_set), biotype, 
                                                                                       gencode)
    barplot(results, out_path, file_name, title_string, adjust_y=False)


def get_highest_cov_alns(cur, genomes):
    """
    Dictionary mapping each genome to a dictionary reporting each highest coverage alignment and its metrics
    """
    return {genome: sql_lib.highest_cov_aln(cur, genome) for genome in genomes}


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.comparativeAnnotationDir)
    highest_cov_dict = get_highest_cov_alns(cur, args.genomes)
    # we need all IDs for this Gencode set to have a consistent denominator
    gencode_ids = seq_lib.get_gp_ids(args.annotationGp)
    # genome_order = find_genome_order(highest_cov_dict, gencode_ids)
    genome_order = hard_coded_genome_order
    for biotype in get_all_biotypes(args.attributePath):
        chr_y_ids = sql_lib.get_ids_by_chromosome(cur)
        biotype_ids = get_all_ids(args.attributePath, biotype=biotype)
        filter_set = (biotype_ids & gencode_ids) - chr_y_ids
        if len(filter_set) > 200:  # hardcoded cutoff to avoid issues where this biotype/gencode mix is nearly empty
            base_file_name = args.gencode
            out_path = os.path.join(args.outDir, biotype)
            mkdir_p(out_path)
            cov_ident_wrapper(highest_cov_dict, genome_order, out_path, base_file_name, biotype, args.gencode, 
                              filter_set)
            cat_plot_wrapper(cur, highest_cov_dict, genome_order, out_path, base_file_name, biotype, args.gencode, 
                             filter_set)
            paralogy_plot(cur, genome_order, out_path, base_file_name, biotype, args.gencode, filter_set)
            num_good_pass(highest_cov_dict, cur, genome_order, out_path, base_file_name, biotype, args.gencode,
                          filter_set)


if __name__ == "__main__":
    main()