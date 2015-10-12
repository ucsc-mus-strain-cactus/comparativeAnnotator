import os

import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib

transcripts = seq_lib.get_gene_pred_transcripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/C57B6NJ/all/transMapGencodeBasicVM4.gp")
transcript_dict = seq_lib.transcript_list_to_dict(transcripts)
annotations = seq_lib.get_gene_pred_transcripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp")
annotation_dict = seq_lib.transcript_list_to_dict(annotations)
alignments = psl_lib.read_psl("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/C57B6NJ/all/transMapGencodeBasicVM4.psl")
alignment_dict = psl_lib.get_psl_dict(alignments)
seq_dict = seq_lib.get_sequence_dict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1509/C57B6NJ.fa")
ref_seq_dict = seq_lib.get_sequence_dict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1509/C57B6J.fa")



augustusTranscripts = seq_lib.get_gene_pred_transcripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/augustus/tmr/C57B6NJ.gp")
augustusTranscriptDict = seq_lib.transcript_list_to_dict(augustusTranscripts)


def compare_intron_to_reference(intron, a, aln, compare_dict, ref_dict):
    """
    For use with the splicing classifiers. Given an index of an intron in t that has a problem,
    determines if this intron exists in the reference. Then, determines if this reference splice
    also has the problem defined by compare_dict. Returns True if the reference also has a splicing problem.

    TODO: I don't understand why I am effectively adding 2 below, once before and once after. a_start cannot ever have
    been None because this would raise a TypeError. I am scared to change this without a test...
    """
    a_start = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.start - 1)) + 1
    a_stop = a.transcript_coordinate_to_chromosome(aln.target_coordinate_to_query(intron.stop))
    for a_intron in a.intron_intervals:
        if a_intron.start == a_start and a_intron.stop == a_stop:
            ref_seq = a_intron.get_sequence(ref_dict, strand=True)
            donor, acceptor = ref_seq[:2], ref_seq[-2:]
            if donor not in compare_dict or compare_dict[donor] != acceptor:
                return True
    return False


canonical = {"GT": "AG"}
for aln_id, aln in alignment_dict.iteritems():
    t = transcript_dict[aln_id]
    a = annotation_dict[psl_lib.remove_alignment_number(aln_id)]
    for intron in t.intron_intervals:
        if len(intron) <= 30:
            continue
        elif not (intron.start >= t.thick_start and intron.stop < t.thick_stop):
            continue
        seq = intron.get_sequence(seq_dict, strand=True)
        donor, acceptor = seq[:2], seq[-2:]
        if donor not in canonical or canonical[donor] != acceptor:
            # is this a intron that exists in the reference that also has this problem?
            if compare_intron_to_reference(intron, a, aln, canonical, ref_seq_dict) is True:
                break


def codon_pair_iterator(a, t, aln, target_seq_dict, query_seq_dict):
    """
    Inputs:
    Transcript objects representing the annotation (query) transcript and the target transcript.
    PslRow object that represents the alignment between the transcript objects.
    SeqDicts/TwoBitFileObjs that contain the genomic sequence for these two transcripts

    Order is (target_cds_pos, target, query)
    """
    target_cds = t.get_cds(target_seq_dict)
    query_cds = a.get_cds(query_seq_dict)
    a_frames = [x for x in a.exon_frames if x != -1]
    if a.strand is True:
        a_offset = a_frames[0]
    else:
        a_offset = 3 - a_frames[-1]
    for i in xrange(a_offset, a.get_cds_length(), 3):
        target_cds_positions = [t.chromosome_coordinate_to_cds(
                                aln.query_coordinate_to_target(
                                a.cds_coordinate_to_transcript(j)))
                                for j in xrange(i, i + 3)]
        if None in target_cds_positions:
            continue
        target_codon = target_cds[target_cds_positions[0]:target_cds_positions[0] + 3]
        query_codon = query_cds[i:i + 3]
        yield target_cds_positions[0], target_codon, query_codon



for aln_id, t in transcript_dict.iteritems():
    a = annotation_dict[psl_lib.remove_alignment_number(aln_id)]
    aln = alignment_dict[aln_id]
    if a.get_cds_length() <= 75 or t.get_cds_length() <= 75:
        continue
    # TODO: this will miss an inframe stop if it is the last 3 bases that are not the annotated stop.
    # use the logic from EndStop to flag this
    codons = list(codon_pair_iterator(a, t, aln, seq_dict, ref_seq_dict))[:-1]