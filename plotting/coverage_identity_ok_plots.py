import argparse
from collections import Counter
import lib.sql_lib as sql_lib
from lib.general_lib import mkdir_p
from plotting.plot_functions import *
from etc.config import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, nargs="+", required=True, help="genomes in this comparison")
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--annotationGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--gencode", type=str, required=True, help="current gencode set being analyzed")
    parser.add_argument("--attributePath", type=str, required=True, help="attribute tsv file")
    return parser.parse_args()


# Hard coded bins used for plots.
paralogy_bins = [0, 1, 2, 3, 4, float('inf')]
identity_bins = [0, 0.0001, 0.995, 0.998, 0.99999999, 1.0]
coverage_bins = [0, 0.0001, 0.8, 0.95, 0.99999999, 1.0]


def paralogy(cur, genome):
    """
    Finds the number of paralogous alignments. This is defined as the number of gene IDs with more than one
    alignment.
    """
    cmd = """SELECT TranscriptId FROM attributes.'{0}'""".format(genome)
    return Counter([x[0] for x in cur.execute(cmd).fetchall()])


def number_categorized(cur, genome, query_fn, biotype=None):
    """
    Finds the alignment IDs categorized by a categorizing function. Can be restricted by biotype
    """
    query = query_fn(genome)
    ok_ids = sql_lib.get_ok_ids(cur, query)
    if biotype is not None:
        biotype_ids = sql_lib.get_biotype_ids(cur, genome, biotype)
        return ok_ids & biotype_ids
    return ok_ids


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


def make_hist(vals, bins, total, g, reverse=False):
    raw = np.histogram(vals, bins)[0]
    if reverse is True:
        raw = raw[::-1]
    norm = raw / (0.01 * total)
    return g, norm, raw


def paralogy_plot(cur, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    results = []
    file_name = "{}_{}".format(base_file_name, "paralogy")
    for g in genome_order:
        p = paralogy(cur, g)
        p = [p.get(x, 0) for x in filter_set]
        g, norm, raw = make_hist(p, paralogy_bins, len(filter_set), g)
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts in {}\nthat have multiple alignments".format(len(filter_set), 
                                                                                                   biotype, gencode)
    legend_labels = ["= {}".format(x) for x in paralogy_bins[:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] 
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
    for query_fn in [alignmentErrors, assemblyErrors]:
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


def num_ok(highest_cov_dict, cur, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    file_name = "{}_numOK".format(base_file_name)
    if biotype == 'protein_coding':
        classifiers = tm_coding_classifiers
    else:
        classifiers = tm_noncoding_classifiers
    results = []
    for genome in genome_order:
        tm_ok = transmap_ok(cur, genome)
        best_ids = set(zip(*highest_cov_dict[genome].itervalues())[0])
        raw = len({x for x in tm_ok if strip_alignment_numbers(x) in filter_set and x in best_ids})
        norm = raw / (0.01 * len(filter_set))
        results.append([genome, norm, raw])
    title_string = "Proportion of {:,} {} transcripts in {}\ncategorized as OK".format(len(filter_set), biotype, 
                                                                                       gencode)
    barplot(results, out_path, file_name, title_string, adjust_y=False)


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.comparativeAnnotationDir)
    highest_cov_dict = {}
    for genome in args.genomes:
        highest_cov_dict[genome] = highest_cov_aln(cur, genome)
    gencode_ids = get_gp_ids(args.annotationGp)
    chr_y_ids = gp_chrom_filter(args.annotationGp)
    # genome_order = find_genome_order(highest_cov_dict, gencode_ids)
    genome_order = hard_coded_genome_order
    for biotype in get_all_biotypes(args.attributePath):
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
            num_ok(highest_cov_dict, cur, genome_order, out_path, base_file_name, biotype, args.gencode, filter_set)


if __name__ == "__main__":
    main()