import os
import lib.sql_lib as sql_lib
import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib


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

