# analyzing the classifiers failed by different subsets of the Gencode annotation
from scripts.consensus import *
from scripts.coverage_identity_ok_plots import *
from scripts.clustering import *
from sonLib.bioio import system
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outDir", required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", required=True, help="directory containing databases")
    parser.add_argument("--basicGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--compGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--attributePath", type=str, required=True, help="attribute tsv file")
    return parser.parse_args()


def get_ids(comp_gp, basic_gp, attr_path):
    comp_ids = get_gp_ids(comp_gp)
    basic_ids = get_gp_ids(basic_gp)
    coding_ids = get_all_ids(attr_path, biotype="protein_coding")
    basic_coding = basic_ids & coding_ids
    comp_coding = comp_ids & coding_ids
    complement_coding = comp_coding - basic_coding
    return basic_coding, comp_coding, complement_coding


def load_data(comp_ann_dir):
    con, cur = attach_databases(comp_ann_dir)
    data = pd.read_sql("select * from C57B6J", con, index_col="ensId")
    data.fillna(0, inplace=True)
    data = data.astype(bool)
    return data


def munge_data(data, filter_set):
    m = data.ix[filter_set]
    m = m[m.sum(axis=1) > 0]  # filter out OK transcripts
    s = m.sum(axis=0)
    normed_s = s / (0.01 * len(m))
    normed_s.sort(ascending=False)
    s.sort(ascending=False)
    drop_low_sums(m, normed_s, cutoff=0.1)
    s = [[x, normed_s[x], y]  for x, y in s.iteritems()]
    return m, s


base_barplot_title = ("Proportion of protein coding transcripts that fail classifiers in the reference\n"
                      "{:,} ({:0.2f}%) not OK transcripts \n")
base_cluster_title = "Hierarchical_clustering_of_classifiers"
file_name_dict = {"Comp": "ref_protein_coding_comprehensive", "Basic": "ref_protein_coding_basic", 
                  "Complement": "ref_protein_coding_complement"}
gencode_dict = {"Comp": "GencodeCompVM4", "Basic": "GencodeBasicVM4", "Complement": "GencodeCompVM4_Specific"}

def main():
    args = parse_args()
    out_path = os.path.join(args.outDir, "reference_classifier_clustering")
    mkdir_p(out_path)
    basic_coding, comp_coding, complement_coding = get_ids(args.compGp, args.basicGp, args.attributePath)
    data = load_data(args.comparativeAnnotationDir)
    for cat, id_set in [["Basic", basic_coding], ["Comp", comp_coding], ["Complement", complement_coding]]:
        m, s = munge_data(data, id_set)
        percent_not_ok = round(100.0 * len(m) / len(id_set), 2)
        barplot_title = base_barplot_title.format(len(m), percent_not_ok)
        file_name = file_name_dict[cat]
        barplot(s, out_path, file_name, barplot_title)
        tmp_csv = "{}.txt".format(cat)
        m.to_csv(tmp_csv)
        genome = "C57B6J"
        biotype = "protein_coding"
        out_cluster_file = os.path.join(out_path, file_name + "_clustering")
        gencode = gencode_dict[cat]
        system("Rscript {}/scripts/cluster.R {} {} {} {} {} {} {} {}".format(os.getcwd(), tmp_csv, base_cluster_title, 
                                                                             genome, len(m), percent_not_ok,  
                                                                             gencode, biotype, out_cluster_file))
        os.remove(tmp_csv)


if __name__ == "__main__":
    main()