os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
import cPickle as pickle
from collections import defaultdict, OrderedDict
from comparativeAnnotator.database_queries import get_row_dict, get_fail_pass_excel_ids, augustus_eval
from pycbio.sys.dataOps import merge_dicts
from pycbio.sys.mathOps import format_ratio
from pycbio.bio.transcripts import GenePredTranscript
from pycbio.bio.intervals import ChromosomeInterval
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers, remove_augustus_alignment_number, \
    aln_id_is_augustus, aln_id_is_transmap
from comparativeAnnotator.database_queries import *
from comparativeAnnotator.database_schema import ref_tables
from comparativeAnnotator.generate_gene_set import *

def build_aln_dict(tx_dict, aug_tx_dict, paralogy_counts):
    """merge different data dicts"""
    r = {}
    for aug_aln_id, aug_t in aug_tx_dict.iteritems():
        t = tx_dict[remove_augustus_alignment_number(aug_aln_id)]
        c = paralogy_counts[aug_aln_id]
        r[aug_aln_id] = (t, aug_t, c)
    return r


args = loadp("v31_args.pickle")
genome = 'SPRET_EiJ'
args.mode = 'transMap'
ref_genome = 'C57B6J'
from pipeline.config import PipelineConfiguration
cfg = PipelineConfiguration(args, args.geneSets[0])
args = cfg.augustus_cfgs[genome]
if mode_is_aug(args.mode):
    gps = load_gps([args.gp, args.augustus_gp])
else:
    gps = load_gps([args.gp])


transcript_gene_map = get_transcript_gene_map(args.query_genome, args.db)
transcript_biotype_map = get_transcript_biotype_map(args.query_genome, args.db)
ref_gps = get_transcript_dict(args.annotation_gp)
biotype_evals = {}
biotype = 'protein_coding'
gp_path = args.geneset_gps[biotype]

gene_transcript_map = get_gene_transcript_map(args.query_genome, args.db, biotype)
ref_gene_sizes = build_gene_sizes(ref_gps, gene_transcript_map, biotype, transcript_biotype_map)
stats = get_db_rows(args.query_genome, args.target_genome, args.db, biotype, args.mode)
excel_ids, pass_specific_ids, fail_ids = get_fail_pass_excel_ids(ref_genome, genome, args.db, biotype,
                                                                     args.filter_chroms)
aug_ids = augustus_eval(ref_genome, genome, args.db, biotype, args.filter_chroms)
id_names = ["fail_ids", "pass_specific_ids", "excel_ids", "aug_ids"]
id_list = [fail_ids, pass_specific_ids, excel_ids, aug_ids]
data_dict = build_data_dict(id_names, id_list, transcript_gene_map, gene_transcript_map)
binned_transcripts = find_best_transcripts(data_dict, stats, args.mode, biotype, gps)


if mode_is_aug(args.mode) and biotype == "protein_coding":
    is_consensus = True
    transcript_evaluation = OrderedDict((x, 0) for x in ["ExcellentTM", "ExcellentAug", "ExcellentTie",
                                                         "PassTM", "PassAug", "PassTie",
                                                         "Fail"])
else:
    is_consensus = False
    transcript_evaluation = OrderedDict((x, 0) for x in ["Excellent", "Pass", "Fail"])


gene_evaluation = OrderedDict((x, 0) for x in ["Excellent", "Pass", "Fail", "NoTransMap"])
longest_rate = OrderedDict((('Longest', OrderedDict((('AddLongest', 0), ('FailAddLongest', 0)))),
                            ('Rescue', OrderedDict((('GeneRescue', 0), ('FailGeneRescue', 0))))))
consensus = []
for gene_id in binned_transcripts:
    ids_included = set()
    categories = set()
    for ens_id in binned_transcripts[gene_id]:
        # evaluate each transcript for a gene
        best_id, category, tie = binned_transcripts[gene_id][ens_id]
        categories.add(category)
        # best_id could be None based on coverage filters
        if category in ["Excellent", "Pass"] and best_id is not None:
            # if a transcript is of high quality, include it
            consensus.append(best_id)
            ids_included.add(best_id)
            s = evaluate_consensus_tx(best_id, category, tie) if is_consensus is True else category
            transcript_evaluation[s] += 1
    # have we included this gene yet? only count if the transcripts included match the gene biotype.
    # this prevents good mappings of retained introns and such being the only transcript.
    biotype_ids_included = {x for x in ids_included if transcript_biotype_map[strip_alignment_numbers(x)] == biotype}
    # there exists Gencode genes where no transcripts have the parent biotype
    gene_in_consensus = True if len(biotype_ids_included) > 0 or \
                                len(ids_included) == len(binned_transcripts[gene_id]) else False
    # if we have, have we included only short transcripts?
    has_only_short_txs = has_only_short(ids_included, ref_gene_sizes[gene_id], gps)
    if has_only_short_txs is True and gene_in_consensus is True:
        # add the single longest transcript for this gene if it passes filters
        longest_id, tie = find_longest_for_gene(binned_transcripts[gene_id], stats, gps, biotype, ids_included)
        if longest_id in consensus:
            continue
        elif longest_id is not None:
            consensus.append(longest_id)
            transcript_evaluation['Fail'] += 1
            longest_rate['Longest']["AddLongest"] += 1
        else:
            longest_rate['Longest']["FailAddLongest"] += 1
    if gene_in_consensus is True:
        # gene is in consensus, evaluate and move on
        s = evaluate_gene(categories)
        gene_evaluation[s] += 1
    else:
        # attempt to add one longest transcript for this failing gene
        longest_id, tie = find_longest_for_gene(binned_transcripts[gene_id], stats, gps, biotype, ids_included)
        if longest_id is None:
            gene_evaluation['NoTransMap'] += 1
            longest_rate['Rescue']["FailGeneRescue"] += 1
        else:
            category = 'Fail'
            consensus.append(longest_id)
            transcript_evaluation[category] += 1
            gene_evaluation[category] += 1
            longest_rate['Rescue']["GeneRescue"] += 1



select AlignmentId,AlignmentCoverage,AlignmentIdentity from CAROLI_EiJ_Attributes where TranscriptId = 'ENSMUST00000165289.7';
select AugustusAlignmentId,AugustusAlignmentCoverage,AugustusAlignmentIdentity from CAROLI_EiJ_AugustusAttributes where TranscriptId = 'ENSMUST00000165289.7';


ref = initialize_session(args.query_genome, args.db, ref_tables)
gps = load_gps([args.gp, args.augustus_gp])
transcript_gene_map = get_transcript_gene_map(args.query_genome, args.db)
ref_gene_intervals = build_ref_intervals(ref, args.query_genome, args.db)
tgt_intervals = build_tgt_intervals(gps)
biotype = 'protein_coding'

excel_ids, pass_specific_ids, fail_ids = get_fail_pass_excel_ids(ref_genome, genome, args.db, biotype=biotype, filter_chroms=["Y", "chrY"], best_cov_only=True)
aug_pass_ids = augustus_eval(ref_genome, genome, args.db, biotype=biotype, filter_chroms=["Y", "chrY"])
aug_stats = get_row_dict(ref_genome, genome, args.db, 'AugustusTMR', biotype)


n = {}
for aug_aln_id, r in aug_stats.iteritems():
    if aug_aln_id in aug_pass_ids:
        tx_id = r.AlignmentId
        if tx_id in excel_ids:
            category = 'Excellent'
        elif tx_id in pass_specific_ids:
            category = 'Pass'
        elif tx_id in fail_ids:
            category = 'Fail'
        else:
            continue
        n[aug_aln_id] = {'AlignmentCoverage': r.AlignmentCoverage, 'AugustusAlignmentCoverage': r.AugustusAlignmentCoverage,
                     'AlignmentIdentity': r.AlignmentIdentity, 'AugustusAlignmentIdentity': r.AugustusAlignmentIdentity,
                     'Category': category}









import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
df = pd.DataFrame.from_dict(n)
df = df.transpose()
df.dropna(inplace=True)
p = sns.pairplot(df, hue='Category', hue_order=["Excellent", "Pass", "Fail"], plot_kws=dict(s=3.5, alpha=0.2))
p.savefig('nj_ident_cov_plots.png', format='png')
plt.close('all')



for cat in ["Excellent", "Pass", "Fail"]:
    df2 = df[df['Category'] == cat]
    p = sns.pairplot(df2, hue='Category', plot_kws=dict(s=4, alpha=0.5))
    p.savefig('nj_ident_cov_plots_{}.png'.format(cat), format='png')
    plt.close('all')


n = {}
for aug_aln_id, r in aug_stats.iteritems():
    if aug_aln_id in aug_pass_ids:
        n[aug_aln_id] = {'AugustusAlignmentCoverage': r.AugustusAlignmentCoverage,
                     'AugustusAlignmentIdentity': r.AugustusAlignmentIdentity}


df2 = pd.DataFrame.from_dict(n)
df2 = df2.transpose()
df3 = df2[df2.AugustusAlignmentIdentity > 5]
df3 = df3[df3.AugustusAlignmentCoverage > 5]
p = sns.pairplot(df3, diag_kind='kde', plot_kws=dict(s=3, alpha=0.5))
p.map_lower(sns.kdeplot, cmap='Blues_d')
p.savefig('nj_aug_stats.png', format='png')
plt.close('all')



import os
import pandas as pd
import cPickle as pickle

from jobTree.scriptTree.target import Target

from pycbio.sys.introspection import classes_in_module
from pycbio.bio.transcripts import get_transcript_dict, convert_strand
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.bio import get_sequence_dict
from pycbio.sys.dataOps import grouper
from pycbio.sys.fileOps import tmpFileGet, ensureDir
from pycbio.sys.sqliteOps import ExclusiveSqlConnection

from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number, remove_augustus_alignment_number, strip_alignment_numbers
import comparativeAnnotator.classifiers as classifiers
import comparativeAnnotator.alignment_classifiers as alignment_classifiers
import comparativeAnnotator.augustus_classifiers as augustus_classifiers
import comparativeAnnotator.alignment_attributes as alignment_attributes

args = loadp("mouse_args.pickle")
genome = 'C57BL_6NJ'
ref_genome = 'C57B6J'
from pipeline.config import PipelineConfiguration
cfg = PipelineConfiguration(args, args.geneSets[0])
args = cfg.query_target_cfgs[genome].comp_ann



def build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_counts, coverage_recs):
    """merge different data dicts"""
    r = {}
    for aln_id, aln in psl_dict.iteritems():
        a = ref_dict[remove_alignment_number(aln_id)]
        t = tx_dict.get(aln_id, None)  # not all alignments have a transcript
        ref_aln = ref_psl_dict[remove_alignment_number(aln_id)]
        c = paralogy_counts[aln_id]
        cov = coverage_recs[aln_id]
        r[aln_id] = (a, t, aln, ref_aln, c, cov)
    return r


ref_dict = get_transcript_dict(args.annotation_gp)
tx_dict = get_transcript_dict(args.target_gp)
psl_dict = get_alignment_dict(args.psl)
ref_psl_dict = get_alignment_dict(args.ref_psl)
seq_dict = get_sequence_dict(args.fasta)
ref_seq_dict = get_sequence_dict(args.ref_fasta)

paralogy_counts = alignment_attributes.paralogy(psl_dict)
coverage_recs = alignment_attributes.highest_cov_aln(psl_dict)
aln_dict = build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_counts, coverage_recs)

ref_fasta = query_seq_dict = get_sequence_dict(args.ref_fasta)
tgt_fasta = ref_seq_dict = get_sequence_dict(args.fasta)
aln_classifier_fns = [x() for x in classes_in_module(alignment_classifiers)]
classifier_fns = [x() for x in  classes_in_module(classifiers)]
r_details = {}


tx_id = 'ENSMUST00000034934.14'
aln_id = 'ENSMUST00000034934.14-1'

recs = []
for aln_id in names_set:
    a = ref_dict[strip_alignment_numbers(aln_id)]
    t = tx_dict[aln_id]
    aln = psl_dict[aln_id]
    if a.thick_start == a.start and find_offset(a.exon_frames, a.strand) == 0 and t.strand is True:
        recs.append([a, t, aln])


cds_filter_fn=lambda intron, t: True
mult3=None
skip_n=False
from comparativeAnnotator.classifiers import *
c = UnknownGap()
c(tx_dict[aln_id], tgt_fasta, cds_filter_fn, mult3, skip_n)


from comparativeAnnotator.database_queries import *
tgt, ref = initialize_tm_session(ref_genome, genome, args.db)
metrics = {}
for aln_id in fail_ids:
    r = tgt_ref_join(tgt, ref)
    r = r.where(tgt.classify.AlignmentId == aln_id)
    q = list(r.naive().execute())[0]
    m = {"NumberIntrons": q.NumberIntrons, "NumberMissingOriginalIntrons": q.NumberMissingOriginalIntrons}
    metrics[aln_id] = m


import pandas as pd
df = pd.DataFrame.from_dict(metrics)
df = df.transpose()
df2 = df[df['NumberIntrons'] <= 10]
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
p = sns.jointplot(x="NumberIntrons", y="NumberMissingOriginalIntrons", data=df, kind="kde")
p.savefig('intron_classifiers.png', format='png')
plt.close('all')


r = tgt_ref_join_aln_id(tgt, ref)
r = intron_inequality(r, tgt, ref)
r = noncoding_classify(r, tgt, ref, coverage=90.0, percent_unknown=5.0)



os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
import argparse
import os
import cPickle as pickle
from collections import defaultdict, OrderedDict
from pycbio.sys.sqliteOps import open_database, get_multi_index_query_dict, get_query_dict
from pycbio.sys.fileOps import ensureFileDir
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map, get_transcript_biotype_map
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.scripts.cgp_consensus import *
genome = 'PWK_PhJ'
ref_genome = 'C57B6J'
plot_args, args_holder = loadp('cgp_args.pickle')
args = args_holder[genome]

transcript_biotype_map = get_transcript_biotype_map(args.ref_genome, args.comp_db)
gene_transcript_map = get_gene_transcript_map(args.ref_genome, args.comp_db, biotype="protein_coding")
transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db, biotype="protein_coding")
con, cur = open_database(args.cgp_db)
# load both consensus and CGP into dictionaries
tmr_consensus_dict = get_transcript_dict(args.consensus_gp)
cgp_dict = get_transcript_dict(args.cgp_gp)
# load the BLAT results from the sqlite database
cgp_stats_query = "SELECT CgpId,EnsId,AlignmentCoverage,AlignmentIdentity FROM '{}_cgp'".format(args.genome)
cgp_stats_dict = get_multi_index_query_dict(cur, cgp_stats_query, num_indices=2)
consensus_stats_query = ("SELECT EnsId,AlignmentCoverage,AlignmentIdentity FROM "
                         "'{}_consensus'".format(args.genome))
consensus_stats_dict = get_query_dict(cur, consensus_stats_query)
# load the intron bits
intron_dict = load_intron_bits(args.cgp_intron_bits)
# final dictionaries
cgp_consensus = {}
metrics = {}
# save all CGP transcripts which have no associated genes
find_new_transcripts(cgp_dict, cgp_consensus, metrics)
# save all CGP transcripts whose associated genes are not in the consensus
consensus_genes = {x.name2 for x in tmr_consensus_dict.itervalues()}
find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, cgp_consensus, metrics, tmr_consensus_dict,
                         gene_transcript_map)
# remove all such transcripts from the cgp dict before we evaluate for updating
cgp_dict = {x: y for x, y in cgp_dict.iteritems() if x not in cgp_consensus}
update_transcripts(cgp_dict, tmr_consensus_dict, transcript_gene_map, intron_dict, cgp_consensus, metrics,
                   cgp_stats_dict, consensus_stats_dict, transcript_biotype_map)
deduplicated_consensus = deduplicate_cgp_consensus(cgp_consensus, metrics)
evaluate_cgp_consensus(deduplicated_consensus, metrics)


from pyfasta import Fasta
aln_id = 'ENSMUST00000130350.1'
ref_transcript_fasta = Fasta(args.align_cds.refTranscriptFasta)
target_genome_fasta = Fasta(args.align_cds.targetGenomeFasta)
gp = consensus_dict[aln_id]
cds = gp.get_cds(target_genome_fasta)


args = args.tmr
from comparativeAnnotator.augustus.run_augustus import *
ens_id = 'ENSMUST00000074051.5-10'
ens_id = 'ENSMUST00000179314.2-1'
args.hints_db = 'mouse_v3.db'

fasta = Fasta(args.fasta)
chrom_sizes = {x.split()[0]: x.split()[1] for x in open(args.chrom_sizes)}

for gp_string in open(args.input_gp):
    assert ens_id not in gp_string


gp = GenePredTranscript(gp_string.rstrip().split("\t"))
chrom = gp.chromosome
start = max(gp.start - args.padding, 0)
stop = min(gp.stop + args.padding, chrom_sizes[chrom])
tm_hint = get_transmap_hints(gp_string, args.tm_2_hints_cmd)

rnaseq_hint = get_rnaseq_hints(args.genome, chrom, start, stop, args.hints_db)
hint = "".join([tm_hint, rnaseq_hint])

seq = fasta[chrom][start:stop]
hint_f, seq_f = write_hint_fasta(hint, seq, chrom, './')
cfg_version, cfg_path = args.cfgs.items()[0]
outf_h = open('test.out', 'w')
run_augustus(hint_f, seq_f, gp.name, start, cfg_version, cfg_path, outf_h, gp, args.augustus_bin)
outf_h.close()


tx_dict = get_transcript_dict(args.gp)
ref_tx_dict = get_transcript_dict(args.annotation_gp)
from pycbio.bio.psl import *
psl_dict = get_alignment_dict(args.psl)
from pycbio.bio.bio import *
seq_dict = get_sequence_dict(args.fasta)
ref_seq_dict = get_sequence_dict(args.ref_fasta)
from comparativeAnnotator.comp_lib.annotation_utils import *
for tx_id, t in tx_dict.iteritems():
    aln = psl_dict[tx_id]
    a = ref_tx_dict[strip_alignment_numbers(tx_id)]
    q = list(codon_pair_iterator(a, t, aln, seq_dict, ref_seq_dict))


target_seq_dict = seq_dict
query_seq_dict = ref_seq_dict
target_cds = t.get_cds(target_seq_dict, in_frame=False)  # OOF creates coordinate problems
query_cds = a.get_cds(query_seq_dict, in_frame=False)
a_frames = [x for x in a.exon_frames if x != -1]
a_offset = find_offset(a_frames, a.strand)
q = []
for i in xrange(a_offset, a.cds_size - a.cds_size % 3, 3):
    target_cds_positions = [t.chromosome_coordinate_to_cds(
                            aln.query_coordinate_to_target(
                            a.cds_coordinate_to_transcript(j)))
                            for j in xrange(i, i + 3)]
    if None in target_cds_positions:
        continue
    # sanity check - should probably remove. But should probably write tests too...
    #assert all([target_cds_positions[2] - target_cds_positions[1] == 1, target_cds_positions[1] -
    #            target_cds_positions[0] == 1, target_cds_positions[2] - target_cds_positions[0] == 2])
    target_codon = target_cds[target_cds_positions[0]:target_cds_positions[0] + 3]
    query_codon = query_cds[i:i + 3]
    assert len(target_codon) == len(query_codon) == 3, a.name
    q.append((target_codon, query_codon))

import comparativeAnnotator.comp_lib.annotation_utils as utils
from pycbio.bio.bio import *
from pycbio.bio.psl import *
from pycbio.bio.transcripts import *
tx_dict = get_transcript_dict('/hive/users/ifiddes/ihategit/pipeline/mouse_work_v4/C57B6J/GencodeCompVM8/transMap/SPRET_EiJ.gp')
psl_dict = get_alignment_dict('/hive/users/ifiddes/ihategit/pipeline/mouse_work_v4/C57B6J/GencodeCompVM8/transMap/SPRET_EiJ.psl')
ref_tx_dict = get_transcript_dict('/hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.gp')
seq_dict = get_sequence_dict('/hive/users/ifiddes/ihategit/pipeline/mouse_work_v4/C57B6J/GencodeCompVM8/genome_files/SPRET_EiJ.fa')
ref_seq_dict = get_sequence_dict('/hive/users/ifiddes/ihategit/pipeline/mouse_work_v4/C57B6J/GencodeCompVM8/genome_files/C57B6J.fa')
mult3 = []
not_mult3 = []
unknown = []
for t in tx_dict.itervalues():
    for intron in t.intron_intervals:
        if len(intron) == 0:
            continue
        if utils.short_intron(intron):
            if len(intron) % 3 == 0 and utils.is_cds(intron, t) is True:
                mult3.append(t)
                break
            elif len(intron) % 3 != 0 and utils.is_cds(intron, t) is True:
                not_mult3.append(t)
                break
            elif "N" in intron.get_sequence(seq_dict):
                unknown.append(t)


for aln_id, t in tx_dict.iteritems():
    a = ref_tx_dict[strip_alignment_numbers(aln_id)]
    aln = psl_dict[aln_id]
    q = list(utils.insertion_iterator(a, aln, mult3=True))


for aln_id, t in tx_dict.iteritems():
    a = ref_tx_dict[strip_alignment_numbers(aln_id)]
    aln = psl_dict[aln_id]
    q = list(utils.deletion_iterator(t, aln))
