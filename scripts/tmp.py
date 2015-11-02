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
#aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
ref_aln_dict = psl_lib.get_alignment_dict(ref_psl)

con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/comparativeAnnotation/2015-10-12/GencodeBasicV23", mode="transMap")

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

for aln_id, aln in aln_dict.iteritems():
    t = tx_dict[aln_id]
    a = ref_dict[psl_lib.strip_alignment_numbers(aln_id)]
    aln_starts_ends = get_adjusted_starts_ends(t, aln)
    count = 0
    for ref_exon in a.exons[1:]:
        r = [aln_start - fuzz_distance <= ref_exon.start <= aln_end + fuzz_distance for aln_start, aln_end in
             aln_starts_ends]
        if not any(r):
            count += 1
    classify_dict[aln_id] = count


num_introns = {aln_id: len(tx_dict[aln_id].intron_intervals) for aln_id in tx_dict}
ref_introns = {aln_id: len(ref_dict[psl_lib.strip_alignment_numbers(aln_id)].intron_intervals) for aln_id in tx_dict}
num_short_introns = {aln_id: len([x for x in t.intron_intervals if comp_ann_lib.short_intron(x)]) for aln_id, t in tx_dict.iteritems()}

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
g = g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=75)
g.savefig("pairplot_100.png")


df3 = df[df['# Original Introns'] <= 20]
g = sns.pairplot(df3)
g = g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=50)
g.savefig("pairplot_20.png")


df4 = df2[df2['# Original Introns'] > 5]
g = sns.pairplot(df4)
g = g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=50)
g.savefig("pairplot_100max_5min.png")



