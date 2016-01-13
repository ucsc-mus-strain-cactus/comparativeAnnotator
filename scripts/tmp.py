import os
import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.comp_ann_lib as comp_ann_lib
import itertools

gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/transMap/2015-10-06/transMap/gorilla/transMapGencodeBasicV23.gp"
ref_gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/transMap/2015-10-06/data/wgEncodeGencodeBasicV23.gp"
aug_gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/augustus/tmr/gorilla.output.gp"
aln_psl = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/transMap/2015-10-06/transMap/gorilla/transMapGencodeBasicV23.psl"
ref_psl =  "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/transMap/2015-10-06/data/wgEncodeGencodeBasicV23.psl"
ref_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.fa"
target_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/gorilla.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
ref_aln_dict = psl_lib.get_alignment_dict(ref_psl)
seq_dict = seq_lib.get_sequence_dict(target_fasta)
ref_seq_dict = seq_lib.get_sequence_dict(ref_fasta)

con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_2/comparativeAnnotation/2015-10-12/GencodeBasicV23", mode="augustus")

genome = 'gorilla'
ref_genome = 'human'
biotype = 'protein_coding'
filter_chroms = ["Y", "chrY"]

stats = merge_stats(cur, 'gorilla')
highest_cov_dict = sql_lib.highest_cov_aln(cur, "gorilla")
highest_cov_ids = set(zip(*highest_cov_dict.itervalues())[0])
biotype_ids = sql_lib.get_biotype_aln_ids(cur, 'gorilla', 'protein_coding')
highest_cov_ids &= biotype_ids
best_stats = {x: y for x, y in stats.iteritems() if psl_lib.remove_augustus_alignment_number(x) in highest_cov_ids}
best_tm = {x: y for x, y in best_stats.iteritems() if x in highest_cov_ids}
best_aug = {x: y for x, y in best_stats.iteritems() if psl_lib.remove_augustus_alignment_number(x) in highest_cov_ids and x not in highest_cov_ids}
r = {"higher_cov": [], "higher_ident": [], "higher_both": [], "worse": []}
for aug_id in best_aug:
    aug_cov, aug_ident = best_aug[aug_id]
    tm_cov, tm_ident = best_tm[psl_lib.remove_augustus_alignment_number(aug_id)]
    if aug_cov > tm_cov and aug_ident > tm_ident:
        r["higher_both"].append(aug_id)
    elif aug_cov > tm_cov:
        r["higher_cov"].append(aug_id)
    elif aug_ident > tm_ident:
        r["higher_ident"].append(aug_id)
    else:
        r["worse"].append(aug_id)


transcript_gene_map = sql_lib.get_transcript_gene_map(cur, ref_genome, biotype=None, filter_chroms=filter_chroms)
gene_transcript_map = sql_lib.get_gene_transcript_map(cur, ref_genome, biotype=biotype, filter_chroms=filter_chroms)
stats = merge_stats(cur, 'gorilla')
fail_ids, good_specific_ids, pass_ids = sql_lib.get_fail_good_pass_ids(cur, ref_genome, genome, biotype)
aug_query = etc.config.augustusEval(genome, ref_genome)
aug_ids = sql_lib.get_query_ids(cur, aug_query)
id_names = ["fail_ids", "good_specific_ids", "pass_ids", "aug_ids"]
id_list = [fail_ids, good_specific_ids, pass_ids, aug_ids]
data_dict = build_data_dict(id_names, id_list, transcript_gene_map, gene_transcript_map)
binned_transcripts = find_best_transcripts(data_dict, stats)
consensus = find_consensus(binned_transcripts, stats)
aug_cons = {x for x in consensus if 'aug' in x}


transcript_evaluation = OrderedDict((x, []) for x in ["PassTM", "PassAug", "PassTie", "GoodTM", "GoodAug", "GoodTie",
                                                     "FailTM", "FailAug", "FailTie", "NoTransMap"])
gene_evaluation = OrderedDict((x, []) for x in ["Pass", "Good", "Fail", "NoTransMap"])
gene_fail_evaluation = OrderedDict((x, []) for x in ["Fail", "NoTransMap"])
for gene_id in binned_transcripts:
    categories = set()
    for ens_id in binned_transcripts[gene_id]:
        best_id, category, tie = binned_transcripts[gene_id][ens_id]
        categories.add(category)
        s = evaluate_transcript(best_id, category, tie)
        transcript_evaluation[s].append(ens_id)
    s = evaluate_gene(categories)
    gene_evaluation[s].append(gene_id)
    if s == "Fail":
        best_for_gene = find_best_for_gene(binned_transcripts[gene_id], stats)
        s = evaluate_best_for_gene(best_for_gene)
        gene_fail_evaluation[s].append(gene_id)

highest_covs = sql_lib.highest_cov_aln(cur, genome)
highest_cov_map = {x: y[0] for x, y in highest_covs.iteritems()}

gene_fail_ids = set(gene_fail_evaluation["Fail"])
transcript_fail_ids = []
for x in gene_fail_ids:
    t_ids = gene_transcript_map[x]
    for t_id in t_ids:
        if t_id in highest_cov_map:
            a_id = highest_cov_map[t_id]
            transcript_fail_ids.append(a_id)

df = pd.read_sql("Select AlignmentId,{} FROM main.'gorilla'".format(",".join(etc.config.all_classifiers)), con, index_col="AlignmentId")
df2 = df.ix[transcript_fail_ids]
gene_ids = {n: transcript_gene_map[psl_lib.strip_alignment_numbers(n)] for n in df2.index}
df3 = df2.copy()
df3["GeneId"] = pd.Series(gene_ids)
df3.to_csv("failed_gene_classifiers.tsv", sep="\t")

no_tm_ids = set(transcript_evaluation["NoTransMap"])
no_tm_genes = {n: transcript_gene_map[psl_lib.strip_alignment_numbers(n)] for n in no_tm_ids}



def augustusEval(genome):
    query = ("SELECT augustus.'gorilla'.AugustusAlignmentId FROM augustus_attributes.'{0}' JOIN main.'{0}' ON "
             "main.'{0}'.AlignmentId = augustus_attributes.'{0}'.AlignmentId JOIN augustus.'{0}' USING "
             "(AugustusAlignmentId) WHERE (AugustusNotSameStart = 0 OR "
             "(main.'{0}'.HasOriginalStart = 1 OR main.'{0}'.StartOutOfFrame = 1)) AND "
             "(AugustusNotSameStop = 0 OR HasOriginalStop = 1) AND "
             "(AugustusExonGain = 0 OR (main.'{0}'.HasOriginalStart = 1 OR main.'{0}'.HasOriginalStop = 1)) AND "
             "(AugustusNotSimilarTerminalExonBoundaries = 0 OR "
             "(main.'{0}'.HasOriginalStart = 1 OR main.'{0}'.HasOriginalStop = 1 OR main.'{0}'.StartOutOfFrame = 1)) "
             "AND AugustusNotSimilarInternalExonBoundaries = 0 AND AugustusNotSameStrand = 0 AND AugustusExonLoss = 0 "
             "AND AugustusParalogy = 0")
    query = query.format(genome)
    return query


query1 = "SELECT * FROM augustus_attributes.'{0}' JOIN main.'{0}' ON main.'{0}'.AlignmentId = augustus_attributes.'{0}'.AlignmentId JOIN augustus.'{0}' USING (AugustusAlignmentId)".format('gorilla')
df1 = pd.read_sql(query1, con, index_col="AugustusAlignmentId")


query2 = "SELECT * FROM attributes.'{0}' JOIN main.'{0}' ON main.'{0}'.AlignmentId = attributes.'{0}'.AlignmentId".format('gorilla')
df2 = pd.read_sql(query2, con, index_col="AlignmentId")


SELECT augustus.'gorilla'.AugustusAlignmentId FROM augustus.'gorilla' JOIN main.'gorilla' ON main.'gorilla'.AlignmentId = augustus.'gorilla'.AlignmentId JOIN main.'human' WHERE (AugustusNotSameStart = 0 OR (main.'gorilla'.HasOriginalStart = 1 OR main.'gorilla'.StartOutOfFrame = 1)) AND (AugustusNotSameStop = 0 OR HasOriginalStop = 1) AND (AugustusExonGain = 0 OR (main.'gorilla'.HasOriginalStart = 1 OR main.'gorilla'.HasOriginalStop = 1)) AND (AugustusNotSimilarTerminalExonBoundaries = 0 OR (main.'gorilla'.HasOriginalStart = 1 OR main.'gorilla'.HasOriginalStop = 1 OR main.'gorilla'.StartOutOfFrame = 1)) AND AugustusNotSimilarInternalExonBoundaries = 0 AND AugustusNotSameStrand = 0 AND AugustusExonLoss = 0 AND AugustusParalogy = 0

("SELECT augustus.'gorilla'.AugustusAlignmentId,augustus_attributes.'gorilla'. FROM augustus.'gorilla' JOIN main.'gorilla' JOIN augustus_stats.'gorilla' "
 "ON main.'gorilla'.AlignmentId = augustus.'gorilla'.AlignmentId WHERE augustus.'gorilla'.AugustusNotSameStart = 1 AND "
 "main.'gorilla'.HasOriginalStart = 0 AND augustus.'gorilla'.AugustusParalogy = 0"

"SELECT TranscriptId FROM main.'human' WHERE BeginStart = 1"

aln_id = "ENST00000308744.10-1"
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]
ref_aln = ref_aln_dict[aln_id[:-2]]
aug_t = aug_dict["augI1-ENST00000308744.10-1"]


aln_id = 'ENST00000340001.8-1'
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]
ref_aln = ref_aln_dict[aln_id[:-2]]



gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/transMap/2015-10-06/transMap/NZO_HlLtJ/transMapGencodeCompVM8.gp"
ref_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/transMap/2015-10-06/data/wgEncodeGencodeCompVM8.gp"
aug_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/augustus/tmr/NZO_HlLtJ.gp"
aln_psl = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/transMap/2015-10-06/transMap/NZO_HlLtJ/transMapGencodeCompVM8.psl"
ref_psl = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/transMap/2015-10-06/data/wgEncodeGencodeCompVM8.psl"
ref_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/assemblies/1509_v2/C57B6J.fa"
target_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/assemblies/1509_v2/NZO_HlLtJ.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
ref_aln_dict = psl_lib.get_alignment_dict(ref_psl)
seq_dict = seq_lib.get_sequence_dict(target_fasta)
ref_seq_dict = seq_lib.get_sequence_dict(ref_fasta)
con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509_v2/comparativeAnnotation/2015-12-17/GencodeCompVM8", mode="augustus")


aln_id = "ENSMUST00000020843.11-1"
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]
ref_aln = ref_aln_dict[aln_id[:-2]]



cmd = "grep {0} /hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/transMap/gorilla/*Comp*psl > /hive/users/ifiddes/example_intron_problems/{0}.psl"
system(cmd.format(aln_id))
cmd2 = "grep {0} /hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/data/*Comp*psl > /hive/users/ifiddes/example_intron_problems/{0}.psl"
system(cmd2.format(psl_lib.strip_alignment_numbers(aln_id)))
system("pslFmt {0}.psl > {0}.psl.txt".format("/hive/users/ifiddes/example_intron_problems/{}".format(aln_id)))
system("pslFmt {0}.psl > {0}.psl.txt".format("/hive/users/ifiddes/example_intron_problems/{}".format(psl_lib.strip_alignment_numbers(aln_id))))
system("grep {0} /hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/transMap/gorilla/*Comp*gp > /hive/users/ifiddes/example_intron_problems/{0}.gp".format(aln_id))
system("genePredFmt /hive/users/ifiddes/example_intron_problems/{0}.gp > /hive/users/ifiddes/example_intron_problems/{0}.gp.txt".format(aln_id))



from collections import defaultdict
classify_dict = {}
details_dict = defaultdict(list)

highest_cov_dict = sql_lib.highest_cov_aln(cur, 'gorilla')
highest_cov_ids = set(zip(*highest_cov_dict.itervalues())[0])


for aln_id, aln in aln_dict.iteritems():
    if aln_id not in highest_cov_ids:
        continue
    t = tx_dict[aln_id]
    a = ref_dict[psl_lib.strip_alignment_numbers(aln_id)]
    aln_starts_ends = comp_ann_lib.get_adjusted_starts_ends(t, aln)
    count = 0
    for ref_exon in a.exons[1:]:
        r = [aln_start - fuzz_distance <= ref_exon.start <= aln_end + fuzz_distance for aln_start, aln_end in
             aln_starts_ends]
        if not any(r):
            count += 1
    classify_dict[aln_id] = count


num_introns = {aln_id: len(tx_dict[aln_id].intron_intervals) for aln_id in highest_cov_ids}
ref_introns = {aln_id: len(ref_dict[psl_lib.strip_alignment_numbers(aln_id)].intron_intervals) for aln_id in highest_cov_ids}
num_short_introns = {aln_id: len([x for x in t.intron_intervals if comp_ann_lib.short_intron(x)]) for aln_id, t in tx_dict.iteritems() if aln_id in highest_cov_ids}

import pandas as pd

df = pd.DataFrame.from_dict({"# Missing Original Introns": classify_dict, "# qStarts": num_introns, "# Original Introns": ref_introns, "# Short Introns": num_short_introns})

df2 = df[df['# Original Introns'] != 0]
df3 = df2[df2['# Missing Original Introns'] != 0]


import cPickle as pickle
with open("df_test.pickle", "w") as outf:
    pickle.dump(df3, outf)


df = pickle.load(open("df_test.pickle"))
df2 = df[df['# Original Introns'] <= 100]
g = sns.pairplot(df2)
g = g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=50)
g.savefig("pairplot_best_cov_100.png")


df3 = df[df['# Original Introns'] <= 20]
g = sns.pairplot(df3)
g = g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=50)
g.savefig("pairplot_best_cov_20.png")



with sns.axes_style("white"):
    g = sns.jointplot(x=df3["# Original Introns"], y=df3["# Missing Original Introns"], kind='kde', color='k')


g.savefig("jointplot_best_cov.png")

import re
r = re.compile("-[0-9]+-")
from collections import Counter
counts = Counter("-".join(r.split(aug_aln_id)) for aug_aln_id in aug_dict if psl_lib.remove_augustus_alignment_number(aug_aln_id) in highest_cov_ids)


file_name = "paralogy_best_aln"
out_path = "./"
paralogy_bins = [0, 1, 2, 3, float('inf')]
norm, raw = make_hist(counts.values(), paralogy_bins, reverse=False, roll=-1)
results = [["intronBits", norm]]
title_string = "AugustusParalogy on best transMaps"
legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                    ["= {}".format(paralogy_bins[0])]
plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


gene_transcript_map = sql_lib.get_gene_transcript_map(cur, 'human', biotype='protein_coding', filter_chroms=['Y', 'chrY'])
transcript_gene_map = sql_lib.get_transcript_gene_map(cur, 'human', biotype='protein_coding', filter_chroms=['Y', 'chrY'])



combined_counts = defaultdict(list)
for aug_aln_id, c in counts.iteritems():
    tx_id = psl_lib.strip_alignment_numbers(aug_aln_id)
    if tx_id not in transcript_gene_map:
        continue
    gene_id = transcript_gene_map[tx_id]
    combined_counts[gene_id].append([aug_aln_id, c])


best_c = {}
for gene_id, vals in combined_counts.iteritems():
    best_c[tx_id] = sorted(vals, key=lambda x: -x[1])[0]



file_name = "paralogy_best_aln_by_gene"
out_path = "./"
paralogy_bins = [1, 2, 3, float('inf')]
values = zip(*best_c.values())[1]
norm, raw = make_hist(values, paralogy_bins, reverse=False, roll=0)
results = [["intronBits", norm]]
title_string = "Fewest number of transcripts produced by Augustus for a given gene"
legend_labels = ["= {}".format(x) for x in paralogy_bins[:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])]
plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)


# filter augustus so we don't have to re-run
save_gp = {}
for gp_string in open("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/augustus/tmr_no_intron_vector/gorilla.output.gp"):
    gp = gp_string.rstrip().split("\t")
    gp = seq_lib.GenePredTranscript(gp)
    tm_gp = tx_dict[psl_lib.remove_augustus_alignment_number(gp.name)]
    if not (gp.start < tm_gp.start or gp.stop > tm_gp.stop):
        save_gp[gp.name] = gp


counts = Counter("-".join(r.split(aug_aln_id)) for aug_aln_id in save_gp if psl_lib.remove_augustus_alignment_number(aug_aln_id) in highest_cov_ids)


save_gp_intron_bits = {}

for gp_string in open("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/augustus/tmr/gorilla.output.gp"):
    gp = gp_string.rstrip().split("\t")
    gp = seq_lib.GenePredTranscript(gp)
    tm_gp = tx_dict[psl_lib.remove_augustus_alignment_number(gp.name)]
    if not (gp.start < tm_gp.start or gp.stop > tm_gp.stop):
        save_gp_intron_bits[gp.name] = gp


bits_counts = Counter("-".join(r.split(aug_aln_id)) for aug_aln_id in save_gp if psl_lib.remove_augustus_alignment_number(aug_aln_id) in highest_cov_ids)


file_name = "paralogy_best_aln_filtered"
out_path = "./"
paralogy_bins = [0, 1, 2, 3, float('inf')]
results = []
for x, y in zip(*[["ShortGaps", "IntronBits"], [counts, bits_counts]]):
    norm, raw = make_hist(y.values(), paralogy_bins, reverse=False, roll=-1)
    results.append([x, norm])
title_string = "AugustusParalogy on best transMaps\noverlap filtered"
legend_labels = ["= {}".format(x) for x in paralogy_bins[1:-2]] + [u"\u2265 {}".format(paralogy_bins[-2])] + \
                    ["= {}".format(paralogy_bins[0])]
plot_lib.stacked_barplot(results, legend_labels, out_path, file_name, title_string)

with open("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/augustus/tmr/gorilla_output_unfiltered.gp") as inf:
    with open("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/augustus/tmr/gorilla.output.gp", "w") as outf:
        for x in inf:
            x = x.rstrip().split("\t")
            if x[0] in save_gp_intron_bits:
                outf.write("\t".join(x) + "\n")


q = "SELECT AlignmentId,main.'{genome}'.UtrUnknownSplice,main.'{genome}'. "



noncoding_classifiers = ["UtrUnknownSplice", "PercentUnknownBases", "HasOriginalIntrons"]
biotype = "snoRNA"
ref_genome = "C57B6J"
for genome in ["C57B6NJ", "PAHARIEiJ"]:
    fail_ids, good_specific_ids, pass_ids = sql_lib.get_fail_good_pass_ids(cur, ref_genome, genome, biotype)
    sql_data = sql_lib.load_data(con, genome, noncoding_classifiers)
    num_original_introns = sql_lib.load_data(con, genome, ["NumberIntrons"], table="attributes")
    query = query.format(genome, biotype, ref_genome)
    df = pd.read_sql(query, con, index_col="AlignmentId")
    df["AlignmentCoverage"] = 100 - df["AlignmentCoverage"]
    munged, stats = munge_intron_data(df, num_original_introns, fail_ids)
    plot_lib.barplot(stats, "./", "{}_barplot".format(genome), "{} Fail {} classifiers".format(genome, biotype))