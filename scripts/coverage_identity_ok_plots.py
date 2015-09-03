import argparse
from scripts.plot_functions import *
from src.queries import assemblyErrors, alignmentErrors


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


def highest_cov_aln(cur, genome):
    """
    Returns the set of alignment IDs that represent the best alignment for each source transcript (that mapped over)
    Best is defined as highest %COV. Also reports the associated coverage value.
    """
    tm_stats = get_tm_stats(cur, genome)  # dictionary mapping each aln_id to [aln_id, %ID, %COV]
    combined_covs = defaultdict(list)
    for aln_id, ident, cov in tm_stats.itervalues():
        tx_id = strip_alignment_numbers(aln_id)
        combined_covs[tx_id].append([aln_id, ident, cov])
    best_cov = {}
    for tx_id, vals in combined_covs.iteritems():
        best_cov[tx_id] = sorted(vals, key=lambda x: -x[2])[0]
    return best_cov


def paralogy(cur, genome):
    """
    Finds the number of paralogous alignments. This is defined as the number of gene IDs with more than one
    alignment.
    """
    cmd = """SELECT TranscriptId FROM attributes.'{0}'""".format(genome)
    return Counter([x[0] for x in cur.execute(cmd).fetchall()])


def number_categorized(cur, genome, cat_fn, biotype=None):
    """
    Finds the alignment IDs categorized by a categorizing function. Can be restricted by biotype
    """
    classify_fields, details_fields, classify_values, classify_operations = cat_fn()
    cmd = """SELECT AlignmentId FROM main.'{0}' JOIN attributes.'{0}' USING (AlignmentId) WHERE (""".format(genome)
    for col, mod in zip(classify_fields[:-1], classify_operations):
        cmd += " {} = ? {}".format(col, mod)
    cmd += " {} = ?)".format(classify_fields[-1])
    if biotype is not None:
        cmd += " AND TranscriptType = '{}'".format(biotype)
    return {x[0] for x in cur.execute(cmd, classify_values).fetchall()}


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
        results.append(make_hist(p, paralogy_bins, len(filter_set), g))
    title_string = "Proportion of {:,} {} transcripts in {}\nthat have multiple alignments".format(len(filter_set), 
                                                                                                   biotype, gencode)
    legend_labels = ["= {}".format(x) for x in paralogy_bins[:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] 
    stacked_barplot(results, legend_labels, out_path, file_name, title_string)


def categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, gencode, filter_set, cat_fn):
    results = []
    for g in genome_order:
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        r = number_categorized(cur, g, cat_fn, biotype=biotype)
        raw = len({x for x in r if strip_alignment_numbers(x) in filter_set and x in best_ids})
        norm = raw / (0.01 * len(filter_set))
        results.append([g, norm, raw])
    title_string = "Proportion of {:,} {} transcripts in {}\ncategorized as {}".format(len(filter_set), biotype, 
                                                                                       gencode, cat_fn.__name__)
    barplot(results, out_path, file_name, title_string, adjust_y=False)


def cat_plot_wrapper(cur, highest_cov_dict, genome_order, out_path, base_file_name, biotype, gencode, filter_set):
    for cat_fn in [assemblyErrors, alignmentErrors]:
        file_name = "{}_{}".format(base_file_name, cat_fn.__name__)
        categorized_plot(cur, highest_cov_dict, genome_order, out_path, file_name, biotype, gencode, filter_set, cat_fn)


def metrics_plot(highest_cov_dict, bins, genome_order, out_path, file_name, biotype, gencode, filter_set, analysis):
    results = []
    for g in genome_order:
        mets = highest_cov_dict[g]
        mets = [eval(analysis) for tx_id, (aln_id, identity, coverage) in mets.iteritems() if tx_id in filter_set]
        mets.extend([0] * (len(filter_set) - len(mets)))
        results.append(make_hist(mets, bins, len(filter_set), g, reverse=True))
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
    for genome in genomes:
        tm_ok = transmap_ok(cur, genome, classifiers)
        best_ids = set(zip(*highest_cov_dict[g].itervalues())[0])
        raw = len({x for x in tm_ok if strip_alignment_numbers(x) in filter_set and x in best_ids})
        norm = raw / (0.01 * len(filter_set))
        results.append([genome, norm, raw])
    title_string = "Proportion of {:,} {} transcripts in {}\ncategorized as OK".format(len(filter_set), biotype, 
                                                                                       gencode)
    barplot(results, out_path, file_name, title_string, adjust_y=False)


def main():
    args = parse_args()
    con, cur = attach_databases(args.comparativeAnnotationDir)
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
            num_ok(highest_cov_dict, cur, genome_order, out_path, base_file_name, biotype, gencode, filter_set)


if __name__ == "__main__":
    main()