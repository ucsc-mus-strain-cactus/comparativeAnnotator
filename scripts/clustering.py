import re
from scripts.consensus import *
from scripts.coverage_identity_ok_plots import *
import pandas as pd
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, getRandomAlphaNumericString

tm_grouped = {"CodingIndels": ["CodingInsertions", "CodingDeletions", "CodingMult3Deletions", "CodingMult3Insertions"],
              "AlignmentGaps": ["CdsGap", "CdsMult3Gap", "UtrGap"],
              "UnknownSplice": ["CdsUnknownSplice", "UtrUnknownSplice"],
              "ContainsN": ["AlignmentAbutsUnknownBases", "UnknownBases", "UnknownGap"]}


tm_noncoding_grouped = {"AlignmentGaps": ["UtrGap"], "UnknownSplice": ["UtrUnknownSplice"],
                        "ContainsN": ["AlignmentAbutsUnknownBases", "UnknownBases", "UnknownGap"]}


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, nargs="+", required=True, help="genomes in this comparison")
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--annotationGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--gencode", type=str, required=True, help="current gencode set being analyzed")
    parser.add_argument("--attributePath", type=str, required=True, help="attribute tsv file")
    return parser


def load_data(con, genome, classifiers):
    base_query = "select AlignmentId, {} from main.".format(", ".join(classifiers)) + "'{}'"
    return pd.read_sql_query(base_query.format(genome), con, index_col="AlignmentId")


def drop_low_sums(m, s, cutoff=0.025):
    """
    Drops any classifiers with below cutoff classifications. Cutoff is a percentage of total.
    """
    for c, v in s.iteritems():
        if v < cutoff:
            m.drop(c, inplace=True, axis=1)


def pre_cluster_munged_data(m, grouping):
    for new_field, old_fields in grouping.iteritems():
        m[new_field] = sum(m[x] for x in old_fields)
        m.drop(old_fields, inplace=True, axis=1)


def munge_data(d, filter_set, pre_cluster=False, coding=False):
    """
    Used to munge input data. Can pre-cluster if you want.
    """
    m = d.ix[filter_set]
    m = m.fillna(0)
    if pre_cluster is True:
        if coding is False:
            pre_cluster_munged_data(m, grouping=tm_noncoding_grouped)
        else:
            pre_cluster_munged_data(m, grouping=tm_grouped)
    m = m.astype(bool)
    s = m.sum(axis=0)
    normed_s = s / (0.01 * len(m))
    normed_s.sort(ascending=False)
    s.sort(ascending=False)
    drop_low_sums(m, normed_s)
    s = [[x, normed_s[x], y]  for x, y in s.iteritems()]
    return m, s


def find_aln_id_set(cur, attr_path, ref_gp_path, genome, biotype, classifiers):
    """
    Finds the set of aln_ids for this combination of reference gencode set, biotype, and fail these classifiers
    (are not OK)
    """
    biotype_names = get_all_ids(attr_path, biotype=biotype)  # load all ens_ids for this biotype
    chr_y_names = gp_chrom_filter(ref_gp_path)
    biotype_set = biotype_names - chr_y_names  # filter out ens_ids that are on chromosome Y
    best_covs = set(zip(*highest_cov_aln(cur, genome).itervalues())[0])  # find which aln_ids are the best
    tm_ok_names = transmap_ok(cur, genome, classifiers)  # find the set of OK with these classifiers
    filter_set = {x for x in best_covs if x not in tm_ok_names and strip_alignment_numbers(x) in biotype_set}
    return filter_set, len(biotype_set)


def main_fn(target, comp_ann_path, attr_path, ref_gp_path, gencode, genome, biotype, base_out_path, method):
    base_clust_title = "Hierarchical_clustering_of_transMap_classifiers"
    barplot_title = ("Proportion of transcripts that fail transMap classifiers\ngenome: {}.\t{:,} ({:0.2f}%) not OK "
                    "transcripts \nGencode set: {}\tBiotype: {}")
    out_path = os.path.join(base_out_path, biotype, "clustering", method, genome)
    con, cur = attach_databases(comp_ann_path)
    if biotype == "protein_coding":
        classifiers = tm_coding_classifiers
        coding = True
    else:
        classifiers = tm_noncoding_classifiers
        coding = False
    sql_data = load_data(con, genome, classifiers)
    filter_set, num_biotype = find_aln_id_set(cur, attr_path, ref_gp_path, genome, biotype, classifiers)
    if num_biotype > 25 and len(filter_set) > 10:
        percent_ok = round(100.0 * len(filter_set) / num_biotype, 2)
        if method == "pre_cluster":
            munged, stats = munge_data(sql_data, filter_set, pre_cluster=True, coding=coding)
        else:
            munged, stats = munge_data(sql_data, filter_set, pre_cluster=False, coding=coding)
        mkdir_p(out_path)
        barplot_title = barplot_title.format(genome, len(filter_set), percent_ok, gencode, biotype)
        out_barplot_file = os.path.join(out_path, "barplot{}_{}.pdf".format(genome, biotype))
        barplot(stats, out_path, out_barplot_file, barplot_title)
        q = getRandomAlphaNumericString(10)
        tmp_path = os.path.join(target.getGlobalTempDir(), "{}_tmp.txt".format(q)) # TODO: make this local
        munged.to_csv(tmp_path)
        out_cluster_file = os.path.join(out_path, "clustering_{}_{}.pdf".format(genome, biotype))
        # TODO: why do we have to use my R? 
        system("Rscript {}/scripts/cluster.R {} {} {} {} {} {} {} {}".format(os.getcwd(), tmp_path, base_clust_title, 
                                                                       genome, len(filter_set), percent_ok, gencode, 
                                                                       biotype, out_cluster_file))


def wrapper(target, comp_ann_path, attr_path, ref_gp_path, gencode, genomes, biotypes, base_out_path):
    for genome in genomes:
        for biotype in biotypes:
            for method in ["full"]:#, "pre_cluster"]:
                target.addChildTargetFn(main_fn, args=(comp_ann_path, attr_path, ref_gp_path, gencode, genome, biotype,
                                                   base_out_path, method))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    #biotypes = get_all_biotypes(args.attributePath)
    biotypes = ["protein_coding", "miRNA", "snoRNA", "snRNA", "lincRNA", "processed_pseudogenes", "unprocessed_pseudogenes", "pseudogenes"]
    job_args = (args.comparativeAnnotationDir, args.attributePath, args.annotationGp, args.gencode, args.genomes,
                biotypes, args.outDir)
    i = Stack(Target.makeTargetFn(wrapper, args=job_args)).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from scripts.clustering import *
    main()
