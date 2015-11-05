import os
import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.comp_ann_lib as comp_ann_lib
import itertools

gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/transMap/gorilla/transMapGencodeBasicV23.gp"
ref_gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/data/wgEncodeGencodeBasicV23.gp"
aug_gp = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/augustus/tmr/gorilla.output.gp"
aln_psl = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/transMap/gorilla/transMapGencodeBasicV23.psl"
ref_psl =  "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/transMap/2015-10-06/data/wgEncodeGencodeBasicV23.psl"
ref_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_1/human.fa"
target_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_1/gorilla.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
ref_aln_dict = psl_lib.get_alignment_dict(ref_psl)

con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/comparativeAnnotation/2015-10-12/GencodeBasicV23", mode="transMap")

genome = 'gorilla'
ref_genome = 'human'
biotype = 'protein_coding'
filter_chroms = ["Y", "chrY"]


aln_id = 'ENST00000340001.8-1'
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]
ref_aln = ref_aln_dict[aln_id[:-2]]



gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/transMap/C57B6NJ/transMapGencodeCompVM7.gp"
ref_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/data/wgEncodeGencodeCompVM7.gp"
aug_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/augustus/tmr/C57B6NJ.output.gp"
aln_psl = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/transMap/C57B6NJ/transMapGencodeCompVM7.psl"
ref_psl = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/data/wgEncodeGencodeCompVM7.psl"
ref_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/assemblies/1509/C57B6J.fa"
target_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/assemblies/1509/C57B6NJ.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
ref_aln_dict = psl_lib.get_alignment_dict(ref_psl)
seq_dict = seq_lib.get_sequence_dict(target_fasta)
con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/comparativeAnnotation/2015-10-12/GencodeCompVM7", mode="augustus")


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