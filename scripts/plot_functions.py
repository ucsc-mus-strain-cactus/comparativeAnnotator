import os
import sys
import itertools
import math
import re
import sqlite3 as sql
import numpy as np
from collections import defaultdict, Counter, OrderedDict
from lib.psl_lib import removeAlignmentNumber, removeAugustusAlignmentNumber
from lib.sqlite_lib import attachDatabase
from lib.general_lib import mkdir_p, DefaultOrderedDict

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.backends.backend_pdf as pltBack
# Below are global variables shared by plotting scripts

# this genetic distance is from Joel. There are functions to derive it, but this makes it consistent across plots.
hard_coded_genome_order = ['C57B6NJ', 'NZOHlLtJ', '129S1', 'FVBNJ', 'NODShiLtJ', 'LPJ', 'AJ', 'AKRJ', 'BALBcJ', 'DBA2J',
                            'C3HHeJ', 'CBAJ', 'WSBEiJ', 'CASTEiJ', 'PWKPhJ', 'SPRETEiJ', 'CAROLIEiJ', 'PAHARIEiJ']

# these classifiers define OK for coding transcripts
tm_coding_classifiers = ["CodingInsertions", "CodingDeletions", "CodingMult3Deletions", "CodingMult3Insertions", 
                         "AlignmentPartialMap", "BadFrame", "BeginStart", "UnknownBases", "AlignmentAbutsUnknownBases",
                         "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                         "EndStop", "InFrameStop", "ShortCds", "StartOutOfFrame", "FrameShift"]

# these classifiers define OK for non-coding transcripts
tm_noncoding_classifiers = ["AlignmentPartialMap", "UtrUnknownSplice", "UtrGap", "UnknownGap", "UnknownBases", 
                            "AlignmentAbutsUnknownBases"]


# used for the plots
width = 9.0
height = 6.0
bar_width = 0.45
# paired_palette has two parallel color spectrums and black as the outgroup color
paired_palette = ["#df65b0", "#dd1c77", "#980043", "#a1dab4", "#41b6c4", "#2c7fb8", "#252525"]
# palette is the seaborn colorbind palette
palette = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56b4e9"]

def skip_header(path):
    """
    The attributes file produced by the pipeline has a header. Skip it. Return a open file handle pointing to line 2.
    """
    f_h = open(path)
    _ = f_h.next()
    return f_h


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    """
    return removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))


def get_all_biotypes(attr_path):
    """
    Returns all biotypes in the attribute database.
    """
    return {x.split()[4] for x in skip_header(attr_path)}


def skip_header(path):
    """
    The attributes file produced by the pipeline has a header. Skip it. Return a open file handle pointing to line 2.
    """
    f_h = open(path)
    _ = f_h.next()
    return f_h


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    """
    return removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))


def get_all_biotypes(attr_path):
    """
    Returns all biotypes in the attribute database.
    """
    return {x.split()[4] for x in skip_header(attr_path)}


def transmap_ok(cur, genome, classify_fields):
    """
    Finds all aIds which are 'OK' based on the classify_fields below
    """
    cmd = """SELECT AlignmentId FROM main.'{0}' WHERE (""".format(genome)
    for col in classify_fields[:-1]:
        cmd += " {} = ? {}".format(col, "AND")
    cmd += " {} = ? )".format(classify_fields[-1])
    vals = [0] * len(classify_fields)
    return {x[0] for x in cur.execute(cmd, vals).fetchall()}


def augustus_ok(cur, genome):
    """
    Finds all aug_aIds which are 'OK' as defined by the fields in classifyFields
    """
    classifyFields = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand', 
                      'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries', 
                      'AugustusNotSimilarInternalExonBoundaries']
    cmd = """SELECT augustus.'{0}'.'AlignmentId' FROM augustus.'{0}' WHERE (""".format(genome)
    for col in classifyFields[:-1]:
        cmd += " augustus.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " augustus.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
    vals = [0] * len(classifyFields)
    return {x[0] for x in cur.execute(cmd, vals).fetchall()}


def get_all_ok(cur, genome, tm_classifiers):
    """
    Adapter function to return the combined sets of ok from augustus and transmap
    """
    return augustus_ok(cur, genome) | transmap_ok(cur, genome, tm_classifiers)


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


def get_all_ids(attr_path, biotype=None, filter_set=set(), id_type="Transcript"):
    """
    returns the set of ensembl IDs in the entire Gencode database pulled from the attribute
    """
    assert id_type in ["Transcript", "Gene"]
    if id_type == "Transcript":
        if biotype is None:
            return {x.split()[3] for x in skip_header(attr_path) if x not in filter_set}
        else:
            return {x.split()[3] for x in skip_header(attr_path) if x.split()[4] == biotype if x not in filter_set}
    else:
        if biotype is None:
            return {x.split()[0] for x in skip_header(attr_path) if x not in filter_set}
        else:
            return {x.split()[0] for x in skip_header(attr_path) if x.split()[4] == biotype if x not in filter_set}        


def get_gp_ids(gp):
    """
    Get all unique gene IDs from a genePred
    """
    return {x.rstrip().split("\t")[0] for x in open(gp)}


def gp_chrom_filter(gp, filter_chrom=re.compile("(Y)|(chrY)")):
    """
    Takes a genePred and lists all transcripts that match filter_chrom
    """
    return {x.rstrip().split("\t")[0] for x in open(gp) if filter_chrom.match(x.rstrip().split("\t")[1])}


def load_gps(gp_paths):
    """
    Get a dictionary mapping all gene IDs from a genePred into its entire record. If the gene IDs are not unique
    this function will not work like you want it to.
    """
    return {l.split()[0]: l for p in gp_paths for l in open(p)}


def get_reverse_name_map(cur, genome, blacklist=set(), whitelist=set(), has_augustus=False):
    """
    creates a dictionary mapping each Gencode ID to all IDs produced by Augustus and transMap
    """
    if len(whitelist) > 0:  # if we have a whitelist, we filter out all blacklist items to make a new whitelist
        whitelist = whitelist - blacklist
    reverse_name_map = defaultdict(list)
    base_cmd = "SELECT {0}.'{1}'.'AlignmentId' FROM {0}.'{1}'"
    aug_cmd = base_cmd.format("augustus", genome)
    tm_cmd = base_cmd.format("main", genome)
    if has_augustus:
        aug_r = cur.execute(aug_cmd).fetchall()
    else:
        aug_r = []
    tm_r = cur.execute(tm_cmd).fetchall()
    for aln_id in itertools.chain(aug_r, tm_r):
        aln_id = aln_id[0]
        ens_id = strip_alignment_numbers(aln_id)
        if (len(whitelist) != 0 and ens_id in whitelist) and ens_id not in blacklist:
            reverse_name_map[ens_id].append(aln_id)
    return reverse_name_map


def transcript_list_to_gene(gene_map, ens_ids):
    """
    Given a set of transcript IDs, returns the matching set of gene IDs
    """
    return {gene_map[strip_alignment_numbers(x)] for x in ens_ids}


def get_gene_biotype_map(attr_path):
    """
    Returns a dictionary mapping all gene IDs to their respective biotypes
    """
    return {x.split()[0]: x.split()[2] for x in skip_header(attr_path)}


def get_gene_map(attr_path):
    """
    Returns a dictionary mapping all transcript IDs to their respective gene IDs
    """
    return {x.split()[3]: x.split()[0] for x in skip_header(attr_path)}


def get_tm_stats(cur, genome):
    """
    Pulls the alignment metrics from the attributes database
    """
    cmd = "SELECT AlignmentId, AlignmentIdentity, AlignmentCoverage FROM attributes.'{}'".format(genome)
    result = cur.execute(cmd).fetchall()
    return {x[0]: x for x in result}


def attach_databases(comp_ann_path, has_augustus=False):
    """
    Attaches all of the databases. Expects comp_ann_path to be the path that comparativeAnnotator wrote to.
    If has_augustus is True, expects this folder to have a augustus database.
    """
    con = sql.connect(os.path.join(comp_ann_path, "classify.db"))
    cur = con.cursor()
    attachDatabase(con, os.path.join(comp_ann_path, "attributes.db"), "attributes")
    if has_augustus:
        assert os.path.exists(os.path.join(comp_ann_path, "augustusClassify.db"))
        attachDatabase(con, os.path.join(comp_ann_path, "augustusClassify.db"), "augustus")
    return con, cur


def make_counts_frequency(counts):
    """
    Convenience function that takes a dict and turns the values into a proportion of the total.
    Returns a list of lists [[name, percent]]. Should probably be a OrderedDict unless you like nonsensical plots.
    """
    normed = OrderedDict()
    tot = sum(counts.values())
    for key, val in counts.iteritems():
        normed[key] = 100.0 * val / tot
    return [[x, y, counts[x]] for x, y in normed.iteritems()]
    return list(normed.iteritems()), list(counts.iteritems())


def init_image(out_folder, comparison_name, width, height):
    pdf = pltBack.PdfPages(os.path.join(out_folder, comparison_name + ".pdf"))
    # width by height in inches
    fig = plt.figure(figsize=(width, height), dpi=300, facecolor='w')
    return fig, pdf


def establish_axes(fig, width, height, border=True, has_legend=True):
    """
    Sets up axes. No idea how this works, Dent's code.
    """
    axLeft = 1.1 / width
    if border is True:
        if has_legend is True:
            axRight = 1.0 - (1.5 / width)
        else:
            axRight = 1.0 - (1.15 / width)
    else:
        if has_legend is True:
            axRight = 1.1 - (1.5 / width)
        else:
            axRight = 1.1 - (1.15 / width)
    axWidth = axRight - axLeft
    axBottom = 1.4 / height
    axTop = 0.90 - (0.4 / height)
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


def adjust_x_labels(ax, names, cutoff1=10, cutoff2=17, cutoff3=25):
    """
    If your xaxis labels have a variable amount of text, this can adjust them individually
    """
    for n, t in itertools.izip(*[names, ax.xaxis.get_major_ticks()]):
        if cutoff2 > len(n) > cutoff1:
            t.label1.set_fontsize(8)
        elif cutoff3 > len(n) >= cutoff2:
            t.label1.set_fontsize(7)
        elif len(n) >= cutoff2:
            t.label1.set_fontsize(6)


def base_barplot(max_y_value, names, out_path, file_name, title_string, border=True, has_legend=True):
    """
    Used to initialize either a stacked or unstacked barplot. Expects the max y value to be somewhere in the 10-100
    range or things will get weird.
    """
    assert 10 <= max_y_value <= 100, (max_y_value, names, out_path, file_name, title_string)
    fig, pdf = init_image(out_path, file_name, width, height)
    ax = establish_axes(fig, width, height, border, has_legend)
    plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, max_y_value])
    plt.tick_params(axis='y', labelsize=9)
    plt.tick_params(axis='x', labelsize=9)
    ax.yaxis.set_ticks(np.arange(0.0, int(max_y_value + 1), max_y_value / 10))
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, int(max_y_value + 1), int(max_y_value / 10))])
    ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(names, rotation=60)
    return ax, fig, pdf


def barplot(results, out_path, file_name, title_string, color="#0072b2", border=True, add_labels=True, adjust_y=True):
    """
    Boilerplate code that will produce a unstacked barplot. Expects results to be a list of lists in the form
    [[name1, value1], [name2, value2]]. The values should be normalized between 0 and 100.
    """
    names, values, raw_values = zip(*results)
    if adjust_y is True:
        max_y_value = math.ceil(max(values) / 10.0) * 10
    else:
        max_y_value = 100.0
    ax, fig, pdf = base_barplot(max_y_value, names, out_path, file_name, title_string, border=border, has_legend=False)
    bars = ax.bar(range(len(names)), values, bar_width, color=color)
    if add_labels is True:
        for i, rect in enumerate(bars):
                ax.text(rect.get_x() + bar_width / 2.0, 0.0 + rect.get_height(), raw_values[i], ha='center', 
                        va='bottom', size=6)
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    pdf.close()
    plt.close()


def stacked_barplot(results, legend_labels, out_path, file_name, title_string, color_palette=palette, border=True):
    """
    Boilerplate code that will produce a unstacked barplot. Expects results to be a list of lists of lists in the form
    [[name1, value1], [name2, value2]]. The values should be normalized between 0 and 100. Should be in the same
    order as legend_labels or your legend will be wrong.
    """
    names, values = zip(*results)
    ax, fig, pdf = base_barplot(100.0, names, out_path, file_name, title_string, border=border, has_legend=True)
    bars = []
    cumulative = np.zeros(len(values))
    for i, d in enumerate(np.asarray(values).transpose()):
        bars.append(ax.bar(range(len(values)), d, bar_width, bottom=cumulative, 
                           color=color_palette[i % len(color_palette)],
                           linewidth=0.0, alpha=1.0))
        cumulative += d
    legend = fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=(1,0.8), fontsize=11, 
                        frameon=True, title="Category")
    fig.savefig(pdf, format='pdf')
    pdf.close()
    plt.close()
