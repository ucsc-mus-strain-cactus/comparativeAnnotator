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
ref_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_1/human.fa"
target_fasta = "/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_1/gorilla.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)

con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/comparative/susie_3_1/comparativeAnnotation/2015-10-12/GencodeBasicV23", mode="transMap")

aln_id = 'ENST00000340001.8-1'
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]


gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/transMap/C57B6NJ/transMapGencodeCompVM7.gp"
ref_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/data/wgEncodeGencodeCompVM7.gp"
aug_gp = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/augustus/tmr/C57B6NJ.output.gp"
aln_psl = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/transMap/2015-10-06/transMap/C57B6NJ/transMapGencodeCompVM7.psl"
ref_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/C57B6J.fa"
target_fasta = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/C57B6NJ.fa"

tx_dict = seq_lib.get_transcript_dict(gp)
ref_dict = seq_lib.get_transcript_dict(ref_gp)
aug_dict = seq_lib.get_transcript_dict(aug_gp)
aln_dict = psl_lib.get_alignment_dict(aln_psl)
con, cur = sql_lib.attach_databases("/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/comparativeAnnotation/2015-10-12/GencodeCompVM7", mode="augustus")


aln_id = "ENSMUST00000020843.11-1"
aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[aln_id[:-2]]

result = []
offset = aln.target_coordinate_to_query(t.transcript_coordinate_to_chromosome(0))
target_splices = {x.stop - offset for x in t.exons[:-1]}
query_splices = {x.stop for x in a.exons[:-1]}
not_original_splices = target_splices - query_splices
for exon in t.exons[:-1]:
    target_splice = exon.stop + offset
    if target_splice in not_original_splices:
        result.append(0)
    else:
        result.append(1)


sorted(target_splices)
sorted(query_splices)