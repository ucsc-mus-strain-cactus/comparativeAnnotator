# analyzing the classifiers failed by different subests of the Gencode annotation

import os
import sys
import itertools
import math
import sqlite3 as sql
from lib.psl_lib import removeAlignmentNumber, removeAugustusAlignmentNumber
from lib.sqlite_lib import attachDatabase
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import pandas as pd
import matplotlib.backends.backend_pdf as pltBack
from scripts.coverage_identity_ok_plots import init_image, establish_axes
from scripts.consensus import get_gp_ids, get_all_ids
from sonLib.bioio import system


comp_ann_dir = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/comparativeAnnotation/2015-08-10_Augustus"
basic_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp"
comp_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeCompVM4.gp"
attr_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
table_dump_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-08-10/GencodeCompVM4/simpleChain/ref_errors.csv"


comp_ids = get_gp_ids(comp_gp)
basic_ids = get_gp_ids(basic_gp)
coding_ids = get_all_ids(attr_path, biotype="protein_coding")
basic_coding = basic_ids & coding_ids
comp_coding = comp_ids & coding_ids
complement_coding = comp_coding - basic_coding

data = pd.read_table(table_dump_path, sep=",", header=0, index_col=0, na_values="")
data.fillna(0, inplace=True)
# merging the gap classifiers and dropping those with near-zero counts
data["SmallIntrons"] = data["CdsGap"] + data["UtrGap"]
data.drop(["ScaffoldGap", "UnknownCdsBases", "UnknownBases", "UnknownGap", "AbutsUnknownBases", "CdsGap", "UtrGap"], inplace=True, axis=1)
data = data.astype(bool)
sums = data.sum(axis=0)
sums.sort(ascending=False)

base_title_string = "Proportion of {:,} {}\nClassified as problematic in the reference"
title_string_dict = {"Comp": "Protein-coding {} in GencodeCompVM4", 
                     "Basic": "Protein-coding {} in GencodeBasicVM4", 
                     "Complement": "Protein-coding {}\nin GencodeCompVM4 and NOT in GencodeBasicVM4"}
file_name_dict = {"Comp": "ref_protein_coding_comprehensive", "Basic": "ref_protein_coding_basic", 
                  "Complement": "ref_protein_coding_complement"}
palette = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56b4e9"]

def establish_axes(fig, width, height, border=True, has_legend=True):
    """
    Sets up axes. No idea how this works, Dent's code.
    """
    axLeft = 1.1 / width
    if border is True:
        if has_legend is True:
            axRight = 0.98 - (1.5 / width)
        else:
            axRight = 0.98 - (1.15 / width)
    else:
        if has_legend is True:
            axRight = 1.1 - (1.5 / width)
        else:
            axRight = 1.1 - (1.15 / width)
    axWidth = axRight - axLeft
    axBottom = 1.2 / height
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



r_base_title_string = "Hierarchical clustering of classifiers\n in {} {}"

bar_width = 0.4
for cat, id_set in [["Basic", basic_coding], ["Comp", comp_coding], ["Complement", complement_coding]]:
    this_data = data.ix[list(id_set)]
    sums = this_data.sum(axis=0)
    sums /= len(id_set) * .01
    title_string = base_title_string.format(len(id_set), title_string_dict[cat].format("Transcripts"))
    fig, pdf = init_image(".", file_name_dict[cat], width=8.0, height=4.0)
    ax = establish_axes(fig, width=8.0, height=4.0, border=True, has_legend=False)
    plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    max_val = math.ceil(max(sums) / 10.0) * 10
    ax.set_ylim([0, max_val])
    plt.tick_params(axis='y', labelsize=8)
    plt.tick_params(axis='x', labelsize=8)
    ax.yaxis.set_ticks(np.arange(0.0, int(max_val + 1), max_val / 10))
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, int(max_val + 1), int(max_val / 10))])
    ax.xaxis.set_ticks(np.arange(0, len(sums)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(sums.index, rotation=55)
    bars = ax.bar(range(len(sums)), list(sums), color="#0072b2")
    for i, rect in enumerate(bars):
        ax.text(rect.get_x() + bar_width / 2.0 + 0.1, 0.02 + rect.get_height(), "{:.3f}".format(sums.iloc[i]), ha='center', va='bottom', size=7)
    fig.savefig(pdf, format='pdf')
    pdf.close()
    # lets do some R clustering
    tmp_csv = "{}.txt".format(cat)
    this_data.to_csv(tmp_csv)
    r_title = r_base_title_string.format(len(id_set), title_string_dict[cat].format("Transcripts"))
    r_title = r_title.replace(" ", "_")
    r_title = r_title.replace("\n", "DUMB_PLACEHOLDER")
    os.system("Rscript scripts/cluster.R {} {} {} &".format(tmp_csv, r_title, file_name_dict[cat] + ".clusters.pdf"))