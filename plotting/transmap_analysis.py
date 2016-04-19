import numpy as np
from comparativeAnnotator.database_queries import get_fail_pass_excel_ids, paralogy, get_ref_ids, get_column,\
    get_transcript_gene_map
import pycbio.plotting.plotting as plot_lib
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers

# Hard coded bins used for plots.
paralogy_bins = [0, 1, 2, 3, 4, float('inf')]
coverage_bins = [0, 0.0001, 95.0, 98.0, 99.0, 99.999999, 100.0]
identity_bins = [0, 0.0001, 99.0, 99.5, 99.8, 99.999999, 100.0]


def paralogy_plot(genomes, ref_genome, biotype, path, db_path):
    results = []
    biotype_ids = get_ref_ids(ref_genome, db_path, biotype)
    for g in genomes:
        p = paralogy(g, ref_genome, db_path, biotype)
        # we roll the list backwards one to put 0 on top
        norm, raw = plot_lib.make_hist(p, paralogy_bins, reverse=False, roll=-1)
        results.append([g, norm])
    title_string = "Proportion of {:,} {} transcripts that have multiple alignments"
    title_string = title_string.format(len(biotype_ids), biotype.replace("_", " "))
    legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                    ["= {}".format(paralogy_bins[0])]
    if len(results) > 5:
        plot_lib.stacked_barplot(results, legend_labels, path, title_string)
    else:
        plot_lib.side_by_side_unequal_barplot(results, legend_labels, path, title_string, ylabel='Percentage of transcripts')


def cov_plot(genomes, ref_genome, biotype, path, db_path):
    results = []
    biotype_ids = get_ref_ids(ref_genome, db_path, biotype)
    for genome in genomes:
        r = get_column(genome, ref_genome, db_path, 'tgt.attrs.AlignmentCoverage', biotype, best_only=True)
        r.extend([0] * (len(biotype_ids) - len(r)))  # add no alignment IDs
        norm, raw = plot_lib.make_hist(r, coverage_bins, reverse=True, roll=0)
        results.append([genome, norm])
    title_string = "transMap alignment coverage breakdown for\n{:,} {} transcripts"
    title_string = title_string.format(len(biotype_ids), biotype.replace("_", " "))
    legend_labels = ["= {0:.1f}%".format(coverage_bins[-1])]
    legend_labels.extend(["< {0:.1f}%".format(x) for x in coverage_bins[2:-1][::-1]])
    legend_labels.append("= {0:.1f}%".format(coverage_bins[0]))
    if len(results) > 5:
        plot_lib.stacked_barplot(results, legend_labels, path, title_string)
    else:
        plot_lib.side_by_side_unequal_barplot(results, legend_labels, path, title_string, ylabel='Percentage of transcripts')


def ident_plot(genomes, ref_genome, biotype, path, db_path):
    results = []
    biotype_ids = get_ref_ids(ref_genome, db_path, biotype)
    for genome in genomes:
        r = get_column(genome, ref_genome, db_path, 'tgt.attrs.AlignmentIdentity', biotype, best_only=True)
        r.extend([0] * (len(biotype_ids) - len(r)))  # add no alignment IDs
        norm, raw = plot_lib.make_hist(r, coverage_bins, reverse=True, roll=0)
        results.append([genome, norm])
    title_string = "transMap alignment identity breakdown for\n{:,} {} transcripts"
    title_string = title_string.format(len(biotype_ids), biotype.replace("_", " "))
    legend_labels = ["= {0:.1f}%".format(coverage_bins[-1])]
    legend_labels.extend(["< {0:.1f}%".format(x) for x in coverage_bins[2:-1][::-1]])
    legend_labels.append("= {0:.1f}%".format(coverage_bins[0]))
    if len(results) > 5:
        plot_lib.stacked_barplot(results, legend_labels, path, title_string)
    else:
        plot_lib.side_by_side_unequal_barplot(results, legend_labels, path, title_string, ylabel='Percentage of transcripts')


def num_pass_excel(genomes, ref_genome, biotype, path, db_path, filter_chroms):
    """
    Number of transcripts in each category.
    """
    results = []
    biotype_ids = get_ref_ids(ref_genome, db_path, biotype)
    for genome in genomes:
        excel_ids, pass_specific_ids, fail_ids = get_fail_pass_excel_ids(ref_genome, genome, db_path, biotype=biotype,
                                                                         filter_chroms=filter_chroms)
        num_no_aln = len(biotype_ids) - sum([len(x) for x in [fail_ids, pass_specific_ids, excel_ids]])
        assert num_no_aln >= 0
        raw = np.array([len(excel_ids), len(pass_specific_ids), len(fail_ids), num_no_aln])
        assert all([x >= 0 for x in raw])
        norm = raw / (0.01 * len(biotype_ids))
        results.append([genome, norm])
    title_string = "Proportion of {:,} {} transcripts categorized as Excellent/Pass/Fail"
    title_string = title_string.format(len(biotype_ids), biotype.replace("_", " "))
    legend_labels = ["Excellent", "Pass", "Fail", "NoAln"]
    if len(results) > 5:
        plot_lib.stacked_barplot(results, legend_labels, path, title_string)
    else:
        plot_lib.side_by_side_unequal_barplot(results, legend_labels, path, title_string, ylabel='Percentage of transcripts')


def num_pass_excel_gene_level(genomes, ref_genome, biotype, path, db_path, filter_chroms):
    """
    Number of genes who have at least one transcript in the highest category they have one for.
    """
    results = []
    transcript_gene_map = get_transcript_gene_map(ref_genome, db_path, biotype=biotype)
    num_genes = len(set(transcript_gene_map.values()))
    for genome in genomes:
        excel_ids, pass_specific_ids, fail_ids = get_fail_pass_excel_ids(ref_genome, genome, db_path, biotype=biotype,
                                                                         filter_chroms=filter_chroms)
        excel_genes = {transcript_gene_map[strip_alignment_numbers(x)] for x in excel_ids}
        pass_specific_genes = {transcript_gene_map[strip_alignment_numbers(x)] for x in pass_specific_ids}
        fail_genes = {transcript_gene_map[strip_alignment_numbers(x)] for x in fail_ids}
        num_excel_genes = len(excel_genes)
        num_pass_genes = len(pass_specific_genes - excel_genes)
        num_fail_genes = len(fail_genes - (pass_specific_genes | excel_genes))
        num_no_aln = num_genes - (num_excel_genes + num_pass_genes + num_fail_genes)
        raw = np.array([num_excel_genes, num_pass_genes, num_fail_genes, num_no_aln])
        assert all([x >= 0 for x in raw])
        norm = raw / (0.01 * num_genes)
        results.append([genome, norm])
    title_string = "Proportion of {:,} {} genes \nwith at least one transcript categorized as Excellent/Pass/Fail"
    title_string = title_string.format(num_genes, biotype.replace("_", " "))
    legend_labels = ["Excellent", "Pass", "Fail", "NoAln"]
    if len(results) > 5:
        plot_lib.stacked_barplot(results, legend_labels, path, title_string)
    else:
        plot_lib.side_by_side_unequal_barplot(results, legend_labels, path, title_string, ylabel='Percentage of transcripts')
