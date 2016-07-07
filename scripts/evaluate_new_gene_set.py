import os
import sys
from collections import OrderedDict, Counter
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.plotting.plotting import side_by_side_unequal_barplot, make_hist, stacked_barplot, stacked_side_by_side_unequal_barplot, stacked_unequal_barplot
from pycbio.sys.dataOps import munge_nested_dicts_for_plotting, flatten_list_of_lists
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map

v1_dir = '/hive/users/ifiddes/ihategit/pipeline/mouse_output_v41/combined_gene_sets/for_browser'
v2_dir = '/hive/users/ifiddes/ihategit/pipeline/mouse_output_v5/combined_gene_sets/for_browser'
v1_gps = [os.path.join(v1_dir, x) for x in os.listdir(v1_dir) if x.endswith('.gp')]
v2_gps = [os.path.join(v1_dir, x) for x in os.listdir(v2_dir) if x.endswith('.gp')]
ordered_genomes = 'C57BL_6NJ NZO_HlLtJ NOD_ShiLtJ FVB_NJ LP_J 129S1_SvImJ AKR_J BALB_cJ A_J DBA_2J CBA_J C3H_HeJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ CAROLI_EiJ Pahari_EiJ'.split()
comp_ann_db = './mouse_work_v5/C57B6J/GencodeCompVM8/comparativeAnnotator/classification.db'
transcript_gene_map = get_transcript_gene_map('C57B6J', comp_ann_db, biotype='protein_coding')

# load all transcripts
v1_txs = {os.path.basename(x).split('.')[0]: get_transcript_dict(x) for x in v1_gps}
v2_txs = {os.path.basename(x).split('.')[0]: get_transcript_dict(x) for x in v2_gps}


def calculate_gene_delta(transcript_gene_map, tx_set_1, tx_set_2):
    # we ignore CGP predictions at first
    gene_set_1 = set([transcript_gene_map[x] for x in tx_set_1.iterkeys() if 'jg' not in x and x in transcript_gene_map])
    gene_set_2 = set([transcript_gene_map[x] for x in tx_set_2.iterkeys() if 'jg' not in x and x in transcript_gene_map])
    # bring in jg predictions
    gene_set_1.update([x.split('.')[0] for x in tx_set_1.iterkeys() if 'jg' in x and x in transcript_gene_map])
    gene_set_2.update([x.split('.')[0] for x in tx_set_2.iterkeys() if 'jg' in x and x in transcript_gene_map])
    return gene_set_2 - gene_set_1


new_genes = {}
for g in ordered_genomes:
    genome_v1 = v1_txs[g]
    genome_v2 = v2_txs[g]
    new_genes[g] = sorted(calculate_gene_delta(transcript_gene_map, genome_v1, genome_v2))
    with open(os.path.join('new_genes_v5', g + '.new_genes.txt'), 'w') as outf:
        for l in new_genes[g]:
            outf.write('{}\n'.format(l))




def calculate_tx_per_gene(transcript_gene_map, tx_set, bins):
    tx_size_counter = Counter()
    for tx in tx_set.iterkeys():
        if 'jg' not in tx and tx in transcript_gene_map:
            tx_size_counter[transcript_gene_map[tx]] += 1
    # bring in unmapped genes
    norm, raw = make_hist(tx_size_counter.values(), bins, reverse=False, roll=0)
    return raw


def construct_dict(ordered_genomes):
    d = OrderedDict()
    for g in ordered_genomes:
        d[g] = OrderedDict()
    return d


tx_bins = [0, 2.1, 5.1, float('inf')]
# construct metrics
delta_sizes = construct_dict(ordered_genomes)
transcripts_per_gene_v1 = OrderedDict()
transcripts_per_gene_v2 = OrderedDict()
for genome in ordered_genomes:
    genome_v1 = v1_txs[genome]
    genome_v2 = v2_txs[genome]
    tx_delta = len(genome_v2) - len(genome_v1)
    gene_delta = len(calculate_gene_delta(transcript_gene_map, genome_v1, genome_v2))
    delta_sizes[genome]['Gene'] = gene_delta
    delta_sizes[genome]['Transcript'] = tx_delta
    transcripts_per_gene_v1[genome] = calculate_tx_per_gene(transcript_gene_map, genome_v1, tx_bins)
    transcripts_per_gene_v2[genome] = calculate_tx_per_gene(transcript_gene_map, genome_v2, tx_bins)


# produce plots
results, categories = munge_nested_dicts_for_plotting(delta_sizes)
names, values = zip(*results)
from pycbio.plotting.plotting import *
width = 9.0
height = 6.0
bar_width = 0.45
shorter_bar_width = bar_width / len(values[0])
path = 'delta_sizes_v41_v5.pdf'
border = True
has_legend = True
title_string = 'Change in gene/transcript set sizes between v0.41 and v0.5'
ylabel = u'change in # of genes/transcripts'
fig, pdf = init_image(path, width, height)
ax = establish_axes(fig, width, height, border, has_legend)
plt.tick_params(axis='y', labelsize=9)
plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
plt.tick_params(axis='x', labelsize=9)
ax.set_ylim([-4000, 800])
ax.yaxis.set_ticks([-4000, -3000, -2000, -1000, 0, 200, 400, 600, 800])
ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
ax.xaxis.set_ticklabels(names, rotation=60)
color_palette = palette
bars = []
for i, d in enumerate(np.asarray(values).transpose()):
    bars.append(ax.bar(np.arange(len(values)) + shorter_bar_width * i, d, shorter_bar_width,
                       color=color_palette[i % len(color_palette)], linewidth=0.0, alpha=1.0))


fig.legend([x[0] for x in bars], categories, bbox_to_anchor=(1, 0.8), fontsize=11, frameon=True, title="Category")
if max(len(x) for x in names) > 15:
    adjust_x_labels(ax, names)


fig.savefig(pdf, format='pdf')
plt.close('all')
pdf.close()


legend_labels = [u'\u2264 2', u'\u2264 5', u'> 5']
legend_title = "transcripts/gene"
results = [[x, y] for x, y in transcripts_per_gene_v1.iteritems()]
stacked_unequal_barplot(results, legend_labels, 'transcript_sizes_histogram_v41.pdf', 'Number of transcripts per gene in V0.41',
    ylabel='Number of genes', max_y_value=22000, breaks=8, legend_title=legend_title)


results = [[x, y] for x, y in transcripts_per_gene_v2.iteritems()]
stacked_unequal_barplot(results, legend_labels, 'transcript_sizes_histogram_v5.pdf', 'Number of transcripts per gene in V0.5',
    ylabel='Number of genes', max_y_value=22000, breaks=8, legend_title=legend_title)
