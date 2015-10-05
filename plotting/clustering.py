import re
import argparse
import pandas as pd
from plotting.plot_functions import *
from plotting.coverage_identity_ok_plots import *
from etc.config import *
import lib.sql_lib as sql_lib
from lib.general_lib import mkdir_p
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, getRandomAlphaNumericString


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=str, required=True, help="genome in this comparison")
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--annotationGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--gencode", type=str, required=True, help="current gencode set being analyzed")
    parser.add_argument("--attributePath", type=str, required=True, help="attribute tsv file")
    return parser


def load_data(con, genome, classifiers):
    base_query = "SELECT AlignmentId,{} FROM main.".format(",".join(classifiers)) + "'{}'"
    return pd.read_sql_query(base_query.format(genome), con, index_col="AlignmentId")


def drop_low_sums(m, s, cutoff=2.0):
    """
    Drops any classifiers with below cutoff classifications. Cutoff is a percentage of total.
    """
    for c, v in s.iteritems():
        if v < cutoff:
            m.drop(c, inplace=True, axis=1)


def munge_data(d, filter_set, coding=False):
    """
    Used to munge input data. Can pre-cluster if you want.
    """
    m = d.ix[filter_set]
    m = m.fillna(0)
    m = m.astype(bool)
    s = m.sum(axis=0)
    normed_s = s / (0.01 * len(m))
    normed_s.sort(ascending=False)
    s.sort(ascending=False)
    drop_low_sums(m, normed_s)
    s = [[x, normed_s[x], y]  for x, y in s.iteritems()]
    return m, s


def find_aln_id_set(cur, attr_path, ref_gp_path, genome, biotype, coding):
    """
    Finds the set of aln_ids for this combination of reference gencode set, biotype, and fail these classifiers
    (are not OK)
    """
    biotype_names = get_all_ids(attr_path, biotype=biotype)  # load all ens_ids for this biotype
    chr_y_names = gp_chrom_filter(ref_gp_path)
    biotype_set = biotype_names - chr_y_names  # filter out ens_ids that are on chromosome Y
    best_covs = set(zip(*highest_cov_aln(cur, genome).itervalues())[0])  # find which aln_ids are the best
    tm_ok_names = transmap_ok(cur, genome, coding)  # find the set of OK
    filter_set = {x for x in best_covs if x not in tm_ok_names and strip_alignment_numbers(x) in biotype_set}
    return filter_set, len(biotype_set)


def main_fn(target, comp_ann_path, attr_path, ref_gp_path, gencode, genome, biotype, base_out_path):
    base_clust_title = "Hierarchical_clustering_of_transMap_classifiers"
    base_barplot_title = ("Proportion of transcripts that fail transMap classifiers\ngenome: {}.    {:,} "
                          "({:0.2f}%) not OK transcripts \nGencode set: {}    Biotype: {}")
    out_path = os.path.join(base_out_path, biotype, "clustering", genome)
    con, cur = sql_lib.attach_databases(comp_ann_path)
    if biotype == "protein_coding":
        classifiers = tm_coding_classifiers
        coding = True
    else:
        classifiers = tm_noncoding_classifiers
        coding = False
    sql_data = load_data(con, genome, classifiers)
    filter_set, num_biotype = find_aln_id_set(cur, attr_path, ref_gp_path, genome, biotype, coding)
    if num_biotype > 25 and len(filter_set) > 10:
        percent_not_ok = round(100.0 * len(filter_set) / num_biotype, 2)
        munged, stats = munge_data(sql_data, filter_set, coding=coding)
        mkdir_p(out_path)
        barplot_title = base_barplot_title.format(genome, len(filter_set), percent_not_ok, gencode, biotype)
        out_barplot_file = os.path.join(out_path, "barplot{}_{}".format(genome, biotype))
        barplot(stats, out_path, out_barplot_file, barplot_title)
        # TODO: why can't I use local temp? R fails inexplicably
        tmp_path = os.path.join(target.getGlobalTempDir(), "{}.txt".format(getRandomAlphaNumericString()))
        munged.to_csv(tmp_path)
        out_cluster_file = os.path.join(out_path, "clustering_{}_{}".format(genome, biotype))
        # TODO: why do we have to use my R and export R_HOME?
        system("export R_HOME=/cluster/home/ifiddes/lib64/R && /cluster/home/ifiddes/bin/Rscript {}/plotting/cluster.R "
               "{} {} {} {} {} {} {} {}".format(os.getcwd(), tmp_path, base_clust_title, genome, len(filter_set),
                                                percent_not_ok, gencode, biotype, out_cluster_file))


def wrapper(target, comp_ann_path, attr_path, ref_gp_path, gencode, genome, biotypes, base_out_path):
    for biotype in biotypes:
        target.addChildTargetFn(main_fn, args=(comp_ann_path, attr_path, ref_gp_path, gencode, genome, biotype,
                                               base_out_path))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    biotypes = ["protein_coding", "miRNA", "snoRNA", "snRNA", "lincRNA", "processed_pseudogenes",
                "unprocessed_pseudogenes", "pseudogenes"]
    i = Stack(Target.makeTargetFn(wrapper, args=[args.comparativeAnnotationDir, args.attributePath, args.annotationGp,
                                                 args.gencode, args.genome, biotypes, args.outDir])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from plotting.clustering import *
    main()
