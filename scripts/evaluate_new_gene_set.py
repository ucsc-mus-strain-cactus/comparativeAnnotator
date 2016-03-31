import os
import sys
from collections import OrderedDict, Counter
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.plotting.plotting import side_by_side_unequal_barplot, make_hist, stacked_barplot, stacked_side_by_side_unequal_barplot
from pycbio.sys.dataOps import munge_nested_dicts_for_plotting, flatten_list_of_lists
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map

v1_dir = '/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/comparativeAnnotation/2015-12-17/GencodeCompVM8/cgp_consensus/for_browser'
v2_dir = './mouse_output/CGP_consensus'
v1_gps = [os.path.join(v1_dir, x) for x in os.listdir(v1_dir) if x.endswith('.gp')]
v2_gps = [os.path.join(v2_dir, x) for x in os.listdir(v2_dir) if x.endswith('.gp')]
ordered_genomes = 'C57BL_6NJ NZO_HlLtJ NOD_ShiLtJ FVB_NJ LP_J 129S1_SvImJ AKR_J BALB_cJ A_J DBA_2J CBA_J C3H_HeJ WSB_EiJ CAST_EiJ PWK_PhJ SPRET_EiJ CAROLI_EiJ Pahari_EiJ'.split()
comp_ann_db = './mouse_work/C57B6J/GencodeCompVM8/comparativeAnnotator/classification.db'
transcript_gene_map = get_transcript_gene_map('C57B6J', comp_ann_db, biotype='protein_coding')

# load all transcripts
v1_txs = {os.path.basename(x).split('.')[0]: get_transcript_dict(x) for x in v1_gps}
v2_txs = {os.path.basename(x).split('.')[0]: get_transcript_dict(x) for x in v2_gps}


def calculate_gene_delta(transcript_gene_map, tx_set_1, tx_set_2):
    # we ignore CGP predictions at first
    gene_set_1 = set([transcript_gene_map[x] for x in genome_v1.iterkeys() if 'jg' not in x])
    gene_set_2 = set([transcript_gene_map[x] for x in genome_v2.iterkeys() if 'jg' not in x])
    # bring in jg predictions
    gene_set_1.update([x.split('.')[0] for x in genome_v1.iterkeys() if 'jg' in x])
    gene_set_2.update([x.split('.')[0] for x in genome_v2.iterkeys() if 'jg' in x])
    return len(gene_set_2 - gene_set_1)


def calculate_tx_per_gene(transcript_gene_map, tx_set, bins):
    tx_size_counter = Counter()
    for tx in tx_set.iterkeys():
        if 'jg' not in tx:
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
    gene_delta = calculate_gene_delta(transcript_gene_map, genome_v1, genome_v2)
    delta_sizes[genome] ['Gene']= gene_delta
    delta_sizes[genome]['Transcript'] = tx_delta
    transcripts_per_gene_v1[genome] = calculate_tx_per_gene(transcript_gene_map, genome_v1, tx_bins)
    transcripts_per_gene_v2[genome] = calculate_tx_per_gene(transcript_gene_map, genome_v2, tx_bins)


# produce plots
results, categories = munge_nested_dicts_for_plotting(delta_sizes)
side_by_side_unequal_barplot(results, categories, 'delta_sizes.pdf', 'Change in gene/transcript set sizes between v0.1 and v0.2',
    ylabel=u'change in # of genes/transcripts')

legend_labels = [u'\u2264 2', u'\u2264 5', u'> 5']
legend_title = "transcripts/gene"
results = [[x, y] for x, y in transcripts_per_gene_v1.iteritems()]
stacked_unequal_barplot(results, legend_labels, 'transcript_sizes_histogram_v1.pdf', 'Number of transcripts per gene in V0.1',
    ylabel='Number of genes', max_y_value=22000, breaks=8, legend_title=legend_title)

results = [[x, y] for x, y in transcripts_per_gene_v2.iteritems()]
stacked_unequal_barplot(results, legend_labels, 'transcript_sizes_histogram_v2.pdf', 'Number of transcripts per gene in V0.2',
    ylabel='Number of genes', max_y_value=22000, breaks=8, legend_title=legend_title)
