#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.backends.backend_pdf as pltBack
import numpy as np

import sqlite3 as sql
import os, argparse, re
from collections import defaultdict, Counter, OrderedDict
from itertools import izip

import src.queries
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from lib.general_lib import DirType, functionsInModule, DefaultOrderedDict

color_palette = [( 93, 165, 218),  # m blue
                 (250, 164,  58),  # m orange
                 ( 96, 189, 104),  # m green
                 (241, 124, 167),  # m red
                 (178, 145,  47),  # m brown
                 (178, 118, 178),  # m purple
                 (241,  88,  84),  # m magenta
                 ( 77,  77,  77),  # m grey
                 (222, 207,  63)   # m yellow
                ]  

# put on a 0-1 scale
color_palette = np.asarray(color_palette) / 255.0

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, nargs="+", required=True, help="genomes in this comparison")
    parser.add_argument("--outDir", type=DirType, required=True, help="output directory")
    parser.add_argument("--comparativeAnnotationDir", type=DirType, required=True, help="directory containing databases")
    parser.add_argument("--annotationGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--width", default=8.0, type=float, help="figure width in inches")
    parser.add_argument("--height", default=4.0, type=float, help="figure height in inches")
    parser.add_argument("--biotypes", nargs="+", default=["protein_coding", "lincRNA", "miRNA", "snoRNA", "snRNA", "processed_pseudogene", "pseudogene", "unprocessed_pseudogene"])
    parser.add_argument("--header", type=str, required=True)
    parser.add_argument("--attrs", type=str, required=True, help="attrs")
    args = parser.parse_args()
    return args


def connect_databases(comparativeAnnotationDir):
    con = sql.connect(os.path.join(comparativeAnnotationDir, "classify.db"))
    cur = con.cursor()
    sql_lib.attachDatabase(con, os.path.join(comparativeAnnotationDir, "details.db"), "details")
    sql_lib.attachDatabase(con, os.path.join(comparativeAnnotationDir, "attributes.db"), "attributes")
    return con, cur


def init_image(out_folder, comparison_name, width, height):
    pdf = pltBack.PdfPages(os.path.join(out_folder, comparison_name + ".pdf"))
    # width by height in inches
    fig = plt.figure(figsize=(width, height), dpi=300, facecolor='w')
    return fig, pdf


def get_biotype_transcript_map(attrs, annotation_gp):
    attr_map = {x.split()[3]: x.split()[4] for x in open(attrs)}
    biotype_map = {x.split()[0]: attr_map[x.split()[0]] for x in open(annotation_gp)}
    biotype_counts = Counter(biotype_map.itervalues())
    return biotype_map, biotype_counts


def establish_axes(fig, width, height, border=True):
    """
    Sets up axes. No idea how this works, Dent's code.
    """
    axLeft = 1.1 / width
    if border is True:
        axRight = 0.98 - (1.3 / width)
    else:
        axRight = 1.1 - (1.3 / width)
    axWidth = axRight - axLeft
    axBottom = 0.9 / height
    axTop = 0.9 - (0.4 / height)
    axHeight = axTop - axBottom
    ax = fig.add_axes([axLeft, axBottom, axWidth, axHeight])
    ax.yaxis.set_major_locator(pylab.NullLocator())
    ax.xaxis.set_major_locator(pylab.NullLocator())
    for loc, spine in ax.spines.iteritems():
        if loc in ['left', 'bottom']:
            spine.set_position(('outward', 10))
        elif loc in ['right', 'top']:
            spine.set_color('none')
        else:
            raise ValueError('unknown spine location: %s' % loc)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    return ax


def plot_bars(ax, data, bar_width):
    """
    Data should be a list of lists representing the counts for each bar, (sums to 1)
    """
    bars = []
    cumulative = np.zeros(len(data))
    # transpose data as a matrix to create stacks
    for i, d in enumerate(np.asarray(data).transpose()):
        bars.append(ax.bar(range(len(data)), d, bar_width, bottom=cumulative, color=color_palette[i % len(color_palette)], 
                           linewidth=0.0, alpha=1.0))
        cumulative += d
    return bars


def best_alignment(cur, genome, biotype):
    """
    Returns a set of alignment IDs that have the highest coverage alignments
    """
    b_cmd = """SELECT attributes.'{0}'.AlignmentCoverage, attributes.'{0}'.TranscriptId, attributes.'{0}'.AlignmentId FROM attributes.'{0}' WHERE attributes.'{0}'.GeneType = '{1}'""".format(genome, biotype)
    r = cur.execute(b_cmd).fetchall()
    ids = defaultdict(str)
    cov = defaultdict(float)
    for val, gene, alignment_id in r:
        gene = str(gene)
        if cov[gene] < val:
            cov[gene] = val
            ids[gene] = alignment_id
    return set(ids.values())


def attribute_by_biotype(cur, genome, biotype, attribute):
    """
    Finds the attribute passed for the best alignment per each transcript
    """
    cmd = """SELECT attributes.'{0}'.{2}, attributes.'{0}'.TranscriptId FROM attributes.'{0}' WHERE attributes.'{0}'.'GeneType' = '{1}'""".format(genome, biotype, attribute)
    r = cur.execute(cmd).fetchall()
    t = defaultdict(float)
    for val, gene in r:
        gene = str(gene)
        if t[gene] < val:
            t[gene] = val
    return t


def paralogy(cur, genome):
    """
    Finds the number of paralogous alignments. This is defined as the number of gene IDs with more than one
    alignment.
    """
    cmd = """SELECT attributes.'{0}'.TranscriptId FROM attributes.'{0}'""".format(genome)
    return Counter([x[0] for x in cur.execute(cmd).fetchall()])


def number_categorized(cur, genome, classifyFields, detailsFields, classifyValues, classifyOperations):
    """
    Finds the number of ALIGNMENTS categorized by a categorizing function
    Only done on coding genes right now.
    """
    biotype="protein_coding"
    ids = best_alignment(cur, genome, biotype)
    cmd = """SELECT (AlignmentId) FROM main.'{0}' JOIN attributes.'{0}' USING ('AlignmentId') WHERE (""".format(genome)
    for col, mod in izip(classifyFields[:-1], classifyOperations):
        cmd += " main.'{}'.'{}' = ? {}".format(genome, col, mod)
    cmd += " main.'{0}'.'{1}' = ? AND ATTRIBUTES.'{0}'.'GeneType' = '{2}')".format(genome, classifyFields[-1], biotype)
    if len(detailsFields) == 1:
        cmd += " AND main.'{}'.'{}' = 1".format(genome, detailsFields[0])
    categorized = [x[0] for x in cur.execute(cmd, classifyValues).fetchall() if x[0] in ids]
    return len(categorized), len(ids)


def ok_coding(cur, genome):
    """
    Finds the number of biotype_map whose best alignment has no problems.
    
    OK is defined as biotype_map who do not hit any classifier except for nonsynonymous/synonymous/noncanon splices

    returns the # of OK and the total #
    """
    biotype = "protein_coding"
    classifyFields = ["CodingInsertions","CodingDeletions", "CodingDeletions", "StartOutOfFrame", "FrameShift", "AlignmentAbutsLeft", "AlignmentAbutsRight", 
                      "AlignmentPartialMap", "BadFrame", "BeginStart", "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                      "EndStop", "InFrameStop", "ShortCds", "UnknownBases"]
    cmd = """SELECT main.'{0}'.'AlignmentId' FROM main.'{0}' WHERE (""".format(genome)
    for col in classifyFields[:-1]:
        cmd += " main.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " main.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
    vals = [0] * len(classifyFields)
    # find the best alignment for each transcript
    b_cmd = """SELECT attributes.'{0}'.AlignmentCoverage, attributes.'{0}'.TranscriptId, attributes.'{0}'.AlignmentId FROM attributes.'{0}' WHERE attributes.'{0}'.GeneType = '{1}'""".format(genome, biotype)
    r = cur.execute(b_cmd).fetchall()
    ids = best_alignment(cur, genome, biotype)
    ok_alignments = cur.execute(cmd, vals).fetchall()
    ok_biotype_map = [x[0] for x in cur.execute(cmd, vals).fetchall() if x[0] in ids]
    return len(ok_biotype_map), len(ids)


def plot_stacked_barplot(results, bins, biotype, name, header, out_dir, width, height, num_biotype_map, bar_width=0.4, shorten_name=False):
    # fix attribute name
    if shorten_name == True:
        short_name = re.findall("[A-Z][a-z]+", name)[1].lower()
    else:
        short_name = name
    # make a ratio, reverse order
    results = [(genome, (val / (1.0 * sum(val)))[::-1]) for genome, val in results.iteritems()]
    fig, pdf = init_image(out_dir, header + "_" + biotype + "_" + name, width, height)
    ax = establish_axes(fig, width, height, border=True)
    plt.text(0.5, 1.08, "Proportions of {2} transcripts in biotype '{0}'\nmapped to other strains / species by {1}".format(biotype, short_name, num_biotype_map),
             horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, 1.0])
    plt.tick_params(axis='both', labelsize=8)
    ax.yaxis.set_ticks(np.arange(0.0, 101.0, 10.0) / 100.0)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, 101, 10)])
    ax.xaxis.set_ticks(np.arange(0, len(results)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(zip(*results)[0], rotation=55)
    bars = plot_bars(ax, zip(*results)[1], bar_width)
    bins = bins[::-1]
    legend_labels = ["= {0:.0f}%".format(100.0 * bins[0]), "< {0:.0f}%".format(100.0 * bins[0])] + ["< {}%".format(round(100.0 * x, 3)) for x in bins[2:-2]] + ["= 0%"]
    legend = fig.legend([x[0] for x in bars][::-1], legend_labels[::-1], bbox_to_anchor=(1,0.8), fontsize=11, frameon=True, title=short_name)
    fig.savefig(pdf, format='pdf')
    pdf.close()


def plot_unstacked_barplot(results, out_dir, name, header, width, height, num_biotype_map, bar_width=0.4):
    num_genes = sum(total for num_ok, total in results.itervalues())
    results = [(genome, 1.0 * num_ok / total, num_ok) for genome, (num_ok, total) in results.iteritems()]
    fig, pdf = init_image(out_dir, header + "_" + name, width, height)
    ax = establish_axes(fig, width, height, border=False)
    plt.text(0.5, 1.08, "Proportion of {} successfully transMapped\nprotein coding transcripts that are categorized as {}".format(num_genes, name.split("_")[0]))
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, 1.0])
    plt.tick_params(axis='both', labelsize=8)
    ax.yaxis.set_ticks(np.arange(0.0, 101.0, 10.0) / 100.0)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, 101, 10)])
    ax.xaxis.set_ticks(np.arange(0, len(results)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(zip(*results)[0], rotation=55)
    bars = ax.bar(range(len(results)), zip(*results)[1], bar_width, color=color_palette[0])
    results_num = zip(*results)[2]
    for i, rect in enumerate(bars):
        ax.text(rect.get_x() + bar_width / 2.0, 0.03 + rect.get_height(), results_num[i], ha='center', va='bottom', size=6)
    fig.savefig(pdf, format='pdf')
    pdf.close()    


def paralogy_plot(results, bins, out_dir, header, width, height, bar_width=0.4):
    # make a ratio, reverse order
    results = [(genome, (val / (1.0 * sum(val)))) for genome, val in results.iteritems()]
    fig, pdf = init_image(out_dir, header + "_paralogy", width, height)
    ax = establish_axes(fig, width, height, border=True)
    plt.text(0.5, 1.08, "Number of alignments of protein_coding transcripts",
             horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, 1.0])
    plt.tick_params(axis='both', labelsize=8)
    ax.yaxis.set_ticks(np.arange(0.0, 101.0, 10.0) / 100.0)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, 101, 10)])
    ax.xaxis.set_ticks(np.arange(0, len(results)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(zip(*results)[0], rotation=55)
    bars = plot_bars(ax, zip(*results)[1], bar_width)
    bins = bins[:-1]
    legend_labels = [u"\u2265 {0:.0f}".format(bins[-1] )] + ["= {0:.0f}".format(x) for x in bins[:-1][::-1]]
    legend = fig.legend([x[0] for x in bars], legend_labels[::-1], bbox_to_anchor=(1,0.8), fontsize=10, frameon=True, title="Number of\nAlignments")
    fig.savefig(pdf, format='pdf')
    pdf.close()


def biotype_plot(results, genome, out_dir, header, width, height, bar_width=0.4):
    fig, pdf = init_image(out_dir, header + "_" + genome + "_biotype", width, height)
    ax = establish_axes(fig, width, height, border=True)
    plt.text(0.5, 1.08, "Proportion of reference transcripts with mapping to {} by biotype".format(genome),
             horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    biotypes, mapped, total = zip(*results)
    r = list(np.asarray(map(float, mapped)) / np.asarray(total))
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, 1.0])
    plt.tick_params(axis='both', labelsize=8)
    ax.yaxis.set_ticks(np.arange(0.0, 101.0, 10.0) / 100.0)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, 101, 10)])
    ax.xaxis.set_ticks(np.arange(0, len(results)) + bar_width / 2.0)
    tmp_biotypes = [x.replace("_", "\n") for x in biotypes]
    ax.xaxis.set_ticklabels(tmp_biotypes, rotation=55)
    bars = ax.bar(range(len(results)), r, bar_width, color=color_palette[0])
    for i, rect in enumerate(bars):
        ax.text(rect.get_x() + bar_width / 2.0, 0.03 + rect.get_height(), mapped[i], ha='center', va='bottom', size=6)
    fig.savefig(pdf, format='pdf')
    pdf.close()    


def main():
    args = parse_args()
    con, cur = connect_databases(args.comparativeAnnotationDir)
    identity_bins = [0, 0.0001, 0.995, 0.998, 0.99999999, 1.000001]
    coverage_bins = [0, 0.0001, 0.9, 0.95, 0.99999999, 1.000001]
    paralogy_bins = [1, 2, 3, 4, 100000]
    biotype_map, biotype_counts = get_biotype_transcript_map(args.attrs, args.annotationGp)
    num_coding = biotype_counts['protein_coding']

    # the order of genomes by best average coverage will determine order for all plots and the table
    attribute = "AlignmentCoverage"
    biotype = "protein_coding"
    results = {genome: attribute_by_biotype(cur, genome, biotype, attribute) for genome in args.genomes}
    results = {genome: len([x for x in vals.itervalues() if x > 0.9999999]) for genome, vals in results.iteritems()}
    # sorted by highest attribute to lowest
    genomes = zip(*sorted(results.iteritems(), key = lambda x: -x[1]))[0]

    # stores stats for everything we test to be dumped in a tsv at the end
    statistics = DefaultOrderedDict(list)
    for attribute in ["AlignmentIdentity", "AlignmentCoverage"]:
        if attribute == "AlignmentIdentity":
            bins = identity_bins
        else:
            bins = coverage_bins
        for biotype in args.biotypes:
            results = OrderedDict((genome, attribute_by_biotype(cur, genome, biotype, attribute)) for genome in genomes)
            if any([len(x) == 0 for x in results.itervalues()]):
                continue
            # add in unmapped biotype_map as 0
            for genome, vals in results.iteritems():
                for name in biotype_map:
                    if biotype_map[name] == biotype and name not in vals:
                        vals[name] = 0
            results_hist = OrderedDict((genome, np.histogram(t.values(), bins)[0]) for genome, t in results.iteritems())
            for genome, t in results.iteritems():
                val = [float(y) for x, y in t.iteritems()]
                statistics[biotype + "_" + attribute].append(round(100.0 * sum(val) / len(val), 3))
            num_biotype = biotype_counts[biotype]
            plot_stacked_barplot(results_hist, bins, biotype, attribute, args.header, args.outDir, args.width, 
                                 args.height, num_biotype, shorten_name=True)
    
    for genome in genomes:
        results = OrderedDict((biotype, attribute_by_biotype(cur, genome, biotype, attribute)) for biotype in args.biotypes)
        if any([len(x) == 0 for x in results.itervalues()]):
            continue
        tmp = []
        for biotype, vals in results.iteritems():
            tmp.append([biotype, len([x for x in vals if biotype_map[x] ==  biotype]), biotype_counts[biotype]])
        biotype_plot(tmp, genome, args.outDir, args.header, args.width, args.height)

    results = OrderedDict((genome, ok_coding(cur, genome)) for genome in genomes)
    for genome, (num_ok, total) in results.iteritems():
        statistics["OK_coding"].append(round(100.0 * num_ok / total, 3))
    plot_unstacked_barplot(results, args.outDir, "OK_coding", args.header, args.width, args.height, num_coding)

    results = OrderedDict((genome, paralogy(cur, genome)) for genome in genomes)
    results_hist = OrderedDict((genome, np.histogram(results[genome].values(), paralogy_bins)[0]) for genome in genomes)
    for genome, vals in results_hist.iteritems():
        statistics["paralogy"].append(round(100.0 - (100.0 * vals[0] / sum(vals)), 3))
    paralogy_plot(results_hist, paralogy_bins, args.outDir, args.header, args.width, args.height)

    categories = functionsInModule(src.queries)
    for category in categories:
        detailsFields, classifyFields, classifyValues, classifyOperations = category()
        results = OrderedDict((g, number_categorized(cur, g, classifyFields, detailsFields, classifyValues, 
                                                     classifyOperations)) for g in genomes)
        percent_results = OrderedDict((g, round(100.0 * c / t, 3)) for g, (c, t) in results.iteritems())
        for genome, percent in percent_results.iteritems():
            statistics[category.__name__].append(percent)
        plot_unstacked_barplot(results, args.outDir, category.__name__, args.header, args.width, args.height, num_coding)

    with open(os.path.join(args.outDir, args.outDir, "summary.tsv"), "w") as outf:
        outf.write("genomes\t"+"\t".join(args.genomes)+"\n")
        for category in statistics:
            outf.write(category + "\t" + "\t".join(map(str, statistics[category])) + "\n")        

if __name__ == "__main__":
    main()