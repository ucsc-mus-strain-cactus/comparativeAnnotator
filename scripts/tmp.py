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



def is_fuzzy_intron(intron, aln, ref_starts, offset=5):
    q_gap_start = aln.target_coordinate_to_query(intron.start - 1)
    q_gap_end = aln.target_coordinate_to_query(intron.stop)
    return query_contains_intron(q_gap_start - offset, q_gap_end + offset, ref_starts)


def query_contains_intron(q_gap_start, q_gap_end, ref_starts):
    r = [q_gap_start <= ref_gap <= q_gap_end for ref_gap in ref_starts]
    return True if any(r) else False


def psl_rc(aln):
    aln = copy.deepcopy(aln)
    q_starts = [aln.q_size - (aln.q_starts[i] + aln.block_sizes[i]) for i in xrange(len(aln.q_starts) - 1, -1, -1)]
    t_starts = [aln.t_size - (aln.t_starts[i] + aln.block_sizes[i]) for i in xrange(len(aln.t_starts) - 1, -1, -1)]
    aln.q_starts = q_starts
    aln.t_starts = t_starts
    aln.strand = "+-" if aln.strand == "-" else "-+"
    aln.block_sizes = aln.block_sizes[::-1]
    return aln


[e.start + aln.q_starts[0] for e in t.exons[1:]]


aln = aln_dict[aln_id]
t = tx_dict[aln_id]
a = ref_dict[psl_lib.strip_alignment_numbers(aln_id)]
ref_aln = ref_aln_dict[psl_lib.strip_alignment_numbers(aln_id)]


for aln_id, aln in aln_dict.iteritems():
    t = tx_dict[aln_id]
    ref_aln = ref_aln_dict[psl_lib.strip_alignment_numbers(aln_id)]
    if ref_aln.strand == "-":
        ref_starts = [ref_aln.q_size - (ref_aln.q_starts[i] + ref_aln.block_sizes[i]) for i in
                      xrange(len(ref_aln.q_starts) - 1, -1, -1)]
    else:
        ref_starts = ref_aln.q_starts
    r = [is_fuzzy_intron(intron, aln, ref_starts, 5) for intron in t.intron_intervals]