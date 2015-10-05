import argparse
from plotting.plot_functions import *
from augustus.consensus import 


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--binnedTranscriptDir", required=True)
    parser.add_argument("--attributePath", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--tmGp", required=True)
    parser.add_argument("--compGp", required=True)
    parser.add_argument("--basicGp", required=True)
    return parser.parse_args()


def find_genome_order(binned_transcript_holder, ens_ids):
    """
    Defines a fixed order of genomes based on the most OK protein coding.
    This is deprecated in favor of using a hard coded order provided by Joel.
    """
    ok_counts = []
    for genome, binned_transcripts in binned_transcript_holder.iteritems():
        num_ok = len(binned_transcripts["protein_coding"]["bestOk"])
        ok_counts.append([genome, 1.0 * num_ok / len(ens_ids)])
    genome_order = sorted(ok_counts, key=lambda x: -x[1])
    return zip(*genome_order)[0]


def make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ens_ids, title_string):
    coding_metrics = OrderedDict()
    for g in genome_order:
        bins = binned_transcript_holder[g]['protein_coding']
        coding_metrics[g] = make_counts_frequency(make_tx_counts_dict(bins, filter_set=ens_ids))
    categories = zip(*coding_metrics[g])[0]
    results = [[g, zip(*coding_metrics[g])[1]] for g in coding_metrics]
    stacked_barplot(results, categories, out_path, out_name, title_string, color_palette=paired_palette)


base_title_string = "Proportion of {:,} {}\nOK / not OK In Target Genomes"
title_string_dict = {"Comp": "Protein-coding {} in GencodeCompVM4", 
                     "Basic": "Protein-coding {} in GencodeBasicVM4", 
                     "Complement": "Protein-coding {}\nin GencodeCompVM4 and NOT in GencodeBasicVM4"}
file_name_dict = {"Comp": "protein_coding_comprehensive", "Basic": "protein_coding_basic", 
                  "Complement": "protein_coding_complement"}


def make_coding_transcript_plots(binned_transcript_holder, out_path, comp_gp, basic_gp, attr_path):
    comp_ids = get_gp_ids(comp_gp)
    basic_ids = get_gp_ids(basic_gp)
    coding_ids = get_all_ids(attr_path, biotype="protein_coding")
    basic_coding = basic_ids & coding_ids
    comp_coding = comp_ids & coding_ids
    complement_coding = comp_coding - basic_coding
    #genome_order = find_genome_order(binned_transcript_holder, comp_coding)
    genome_order = hard_coded_genome_order
    for cat, ids in zip(*[["Comp", "Basic", "Complement"], [comp_coding, basic_coding, complement_coding]]):
        title_string = base_title_string.format(len(ids), title_string_dict[cat].format("Transcript"))
        out_name = file_name_dict[cat]
        make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ids, title_string)


def calculate_gene_ok_metrics(bins, gene_map, gene_ids):
    ok_genes = {gene_map[strip_alignment_numbers(x)] for x in bins["bestOk"] if strip_alignment_numbers(x) in gene_map}
    not_ok_genes = {gene_map[strip_alignment_numbers(x)] for x in bins["bestNotOk"] if strip_alignment_numbers(x) in 
                    gene_map and gene_map[strip_alignment_numbers(x)] not in ok_genes}
    fail_genes = gene_ids - (ok_genes | not_ok_genes)
    od = OrderedDict([["Has OK Tx", len(ok_genes)], ["No OK Tx", len(not_ok_genes)], ["No Tx", len(fail_genes)]])
    return make_counts_frequency(od)


def ok_gene_by_biotype(binned_transcript_holder, out_path, attr_path, gene_map, genome_order, biotype):
    biotype_ids = get_all_ids(attr_path, biotype=biotype, id_type="Gene")
    title_string = "Proportion of {:,} {} genes with at least one OK transcript".format(len(biotype_ids), biotype)
    file_name = "{}_gene".format(biotype)
    metrics = OrderedDict()
    for g in genome_order:
        bins = binned_transcript_holder[g][biotype]
        metrics[g] = calculate_gene_ok_metrics(bins, gene_map, biotype_ids)
    categories = zip(*metrics[g])[0]
    results = [[g, zip(*metrics[g])[1]] for g in metrics]
    stacked_barplot(results, categories, out_path, file_name, title_string)


def main():
    plots_path = os.path.join(args.outDir, "geneSetMetrics")
    mkdir_p(plots_path)
    gene_map = get_gene_map(args.attributePath)
    binned_transcript_holder = {}
    for g in os.listdir(args.binnedTranscriptDir):
        binned_transcript_holder[g] = pickle.load(open(os.path.join(args.binnedTranscriptDir, g)))
    make_coding_transcript_plots(binned_transcript_holder, plots_path, args.compGp, args.basicGp, args.attributePath)
    biotype = "protein_coding"
    # ok_gene_by_biotype(binned_transcript_holder, plots_path, args.attributePath, gene_map, args.genome_order, biotype)
    ok_gene_by_biotype(binned_transcript_holder, plots_path, args.attributePath, gene_map, hard_coded_genome_order,
                       biotype)


if __name__ == "__main__":
    main()