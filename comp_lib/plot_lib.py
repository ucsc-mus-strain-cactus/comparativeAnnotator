"""
This file contains convenience functions for plotting.
"""
import os
import itertools
import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.backends.backend_pdf as plt_back
from etc.config import *

__author__ = "Ian Fiddes"

width = 9.0
height = 6.0
bar_width = 0.45


def find_genome_order(highest_cov_dict, biotype_ids):
    """
    Finds the genome order that will be used by all plots. This is deprecated in favor of using a hard coded order
    based on the pipeline config file.
    """
    num_cov = {}
    for g, covs in highest_cov_dict.iteritems():
        num_cov[g] = 1.0 * len({tx_id for tx_id in covs.iterkeys() if tx_id in biotype_ids})
        num_cov[g] /= len(biotype_ids)
    order = sorted(num_cov.iteritems(), key=lambda x: -x[1])
    return zip(*order)[0]


def init_image(out_folder, comparison_name, width, height):
    """
    Sets up a PDF object.
    """
    pdf = plt_back.PdfPages(os.path.join(out_folder, comparison_name + ".pdf"))
    # width by height in inches
    fig = plt.figure(figsize=(width, height), dpi=300, facecolor='w')
    return fig, pdf


def establish_axes(fig, width, height, border=True, has_legend=True):
    """
    Sets up custom axes. The extra space given is adjusted by border and has_legend.
    """
    ax_left = 1.1 / width
    if border is True:
        if has_legend is True:
            ax_right = 1.0 - (1.8 / width)
        else:
            ax_right = 1.0 - (1.15 / width)
    else:
        if has_legend is True:
            ax_right = 1.1 - (1.8 / width)
        else:
            ax_right = 1.1 - (1.15 / width)
    ax_width = ax_right - ax_left
    ax_bottom = 1.4 / height
    ax_top = 0.90 - (0.4 / height)
    ax_height = ax_top - ax_bottom
    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
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


def adjust_x_labels(ax, names, cutoff1=12, cutoff2=18, cutoff3=26):
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


def calculate_y_range(max_y_value, breaks):
    pn = 1.0 * breaks ** math.ceil(math.log10(max_y_value) - 1)
    max_ceil_val = math.ceil(max_y_value / pn) * pn
    return np.arange(0, max_ceil_val + 1, max_ceil_val / breaks)


def base_barplot(max_y_value, names, out_path, file_name, title_string, breaks, border=True, has_legend=True):
    """
    Used to initialize either a stacked or unstacked barplot.
    """
    fig, pdf = init_image(out_path, file_name, width, height)
    ax = establish_axes(fig, width, height, border, has_legend)
    plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, max_y_value])
    plt.tick_params(axis='y', labelsize=9)
    plt.tick_params(axis='x', labelsize=9)
    y_range = calculate_y_range(max_y_value, breaks)
    ax.yaxis.set_ticks(y_range)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in y_range])
    ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(names, rotation=60)
    return ax, fig, pdf


def barplot(results, out_path, file_name, title_string, color="#0072b2", border=True, add_labels=True, adjust_y=True,
            breaks=10.0):
    """
    Boilerplate code that will produce a unstacked barplot. Expects results to be a list of lists in the form
    [[name1, normalized_value1, value1], [name2, normalized_value2, value2]]. Normalized between 0 and 100.
    Assumes that all bars have the same denominator, and so are normalizable.
    """
    names, values, raw_values = zip(*results)
    if adjust_y is True:
        max_y_value = max(values)
    else:
        max_y_value = 100.0
    ax, fig, pdf = base_barplot(max_y_value, names, out_path, file_name, title_string, breaks, border=border, 
                                has_legend=False)
    bars = ax.bar(range(len(names)), values, bar_width, color=color)
    if add_labels is True:
        for i, rect in enumerate(bars):
            v = "{:,}".format(raw_values[i])
            ax.text(rect.get_x() + bar_width / 2.0, 0.0 + rect.get_height(), v, ha='center', va='bottom', size=6)
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()


def stacked_barplot(results, legend_labels, out_path, file_name, title_string, color_palette=palette, border=True,
                    breaks=10.0):
    """
    Boilerplate code that will produce a stacked barplot. Expects results to be a list of lists of lists in the form
    [[name1, [value1a, value1b]], [name2, [value2a, value2b]]. The values should be normalized between 0 and 100.
    Should be in the same order as legend_labels or your legend will be wrong.
    Assumes that all bars have the same denominator, and so are normalizable.
    """
    names, values = zip(*results)
    ax, fig, pdf = base_barplot(100.0, names, out_path, file_name, title_string, breaks, border=border, has_legend=True)
    bars = []
    cumulative = np.zeros(len(values))
    for i, d in enumerate(np.asarray(values).transpose()):
        bars.append(ax.bar(range(len(values)), d, bar_width, bottom=cumulative,
                           color=color_palette[i % len(color_palette)],
                           linewidth=0.0, alpha=1.0))
        cumulative += d
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=(1, 0.8), fontsize=11,
               frameon=True, title="Category")
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()


def base_unequal_barplot(max_y_value, names, out_path, file_name, title_string, ylabel, breaks, border=True,
                         has_legend=True):
    fig, pdf = init_image(out_path, file_name, width, height)
    ax = establish_axes(fig, width, height, border, has_legend)
    plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel(ylabel)
    ax.set_ylim([0, max_y_value])
    plt.tick_params(axis='y', labelsize=9)
    plt.tick_params(axis='x', labelsize=9)
    ax.yaxis.set_ticks(calculate_y_range(max_y_value, breaks))
    ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(names, rotation=60)
    return ax, fig, pdf


def unequal_barplot(results, out_path, file_name, title_string, color_palette=palette, breaks=10.0, border=False, 
                    ylabel="Number of transcripts"):
    """
    Boilerplate code that will produce a barplot. Expects results to be a list of lists of lists in the form
    [[name1, val1], [name2, val2]].
    Should be in the same order as legend_labels or your legend will be wrong.
    """
    names, values = zip(*results)
    max_y_value =  max(values)
    ax, fig, pdf = base_unequal_barplot(max_y_value, names, out_path, file_name, title_string, ylabel, breaks,
                                        border=border, has_legend=True)
    bars = ax.bar(range(len(names)), values, bar_width, color=palette[0])
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()


def stacked_unequal_barplot(results, legend_labels, out_path, file_name, title_string, color_palette=palette,
                            breaks=10.0, border=True, ylabel="Number of transcripts"):
    """
    Boilerplate code that will produce a stacked barplot. Expects results to be a list of lists of lists in the form
    [[name1, [value1a, value1b]], [name2, [value2a, value2b]].
    Should be in the same order as legend_labels or your legend will be wrong.
    """
    names, values = zip(*results)
    max_y_value = max(sum(x) for x in values)
    ax, fig, pdf = base_unequal_barplot(max_y_value, names, out_path, file_name, title_string, ylabel, breaks,
                                        border=border, has_legend=True)
    bars = []
    cumulative = np.zeros(len(values))
    for i, d in enumerate(np.asarray(values).transpose()):
        bars.append(ax.bar(range(len(values)), d, bar_width, bottom=cumulative,
                           color=color_palette[i % len(color_palette)],
                           linewidth=0.0, alpha=1.0))
        cumulative += d
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=(1, 0.8), fontsize=11,
               frameon=True, title="Category")
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()


def side_by_side_unequal_barplot(results, legend_labels, out_path, file_name, title_string, color_palette=palette,
                                 breaks=10.0, border=True, ylabel="Number of transcripts"):
    """
    Boilerplate code that will produce a side by side barplot. Expects results to be a list of lists of lists in the form
    [[name1, [value1a, value1b]], [name2, [value2a, value2b]].
    Should be in the same order as legend_labels or your legend will be wrong.
    """
    names, values = zip(*results)
    shorter_bar_width = bar_width / len(values[0])
    max_y_value = max(sum(x) for x in values)
    ax, fig, pdf = base_unequal_barplot(max_y_value, names, out_path, file_name, title_string, ylabel, breaks,
                                        border=border, has_legend=True)
    bars = []
    for i, d in enumerate(np.asarray(values).transpose()):
        bars.append(ax.bar(np.arange(len(values)) + shorter_bar_width * i, d, shorter_bar_width, 
                           color=color_palette[i % len(color_palette)], linewidth=0.0, alpha=1.0))
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=(1, 0.8), fontsize=11,
               frameon=True, title="Category")
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()


def stacked_side_by_side_unequal_barplot(results, legend_labels, out_path, file_name, title_string, num_columns=2,
                                         color_palette=palette, breaks=10.0, border=True,
                                         ylabel="Number of transcripts"):
    """
    DOESNT WORK YET
    Boilerplate code that will produce a side by side barplot. Expects results to be a list of lists of lists in 
    the form [(colA, colB), (values)]
    At this point only supports pairs, higher orders will not work.
    Should be in the same order as legend_labels or your legend will be wrong.
    """
    def grouper(iterable, n=2):
        return itertools.izip(*[iter(iterable)] * n)
    names, values = zip(*results)
    num_groups = len(names) / num_columns
    num_columns = len(names[0])
    names = [" ".join(x) for x in names]
    shorter_bar_width = bar_width / len(values[0])
    max_y_value = max(sum(x) for x in values)
    ax, fig, pdf = base_unequal_barplot(max_y_value, names, out_path, file_name, title_string, ylabel, breaks,
                                        border=border, has_legend=True)
    rowvals = zip(*values)
    artists = []
    for i, row in enumerate(rowvals):
        cumulative = [np.zeros(num_groups)] * num_columns
        for j in range(num_columns):
            vals = np.array([x for p, x in enumerate(row) if p % num_columns == 0])
            bars = ax.bar(np.arange(len(vals)) + shorter_bar_width * j, vals, shorter_bar_width,
                          color=color_palette[i % len(color_palette)], linewidth=0.0, alpha=1.0)
            cumulative[j] += vals
        artists.append(bars[0])
    fig.legend([x for x in artists[::-1]], legend_labels[::-1], bbox_to_anchor=(1, 0.8), fontsize=11,
               frameon=True, title="Category")
    if max(len(x) for x in names) > 15:
        adjust_x_labels(ax, names)
    fig.savefig(pdf, format='pdf')
    plt.close('all')
    pdf.close()
