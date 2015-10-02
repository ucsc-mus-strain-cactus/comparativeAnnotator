from src.abstract_classifier import Attribute

import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
from lib.general_lib import format_ratio


class TranscriptId(Attribute):
    """
    Creates a column representing the transcript Id
    """
    def run(self):
        self.get_alignment_dict()
        results_dict = {aln_id: psl_lib.remove_alignment_number(aln_id) for aln_id in self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class GeneId(Attribute):
    """
    Creates a column representing the gene Id
    """
    def run(self):
        self.get_attribute_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.attribute_dict[psl_lib.remove_alignment_number(aln_id)].geneId for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class GeneName(Attribute):
    """
    Creates a column representing the gene name
    """
    def run(self):
        self.get_attribute_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.attribute_dict[psl_lib.remove_alignment_number(aln_id)].geneName for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class GeneType(Attribute):
    """
    Creates a column representing the gene type
    """
    def run(self):
        self.get_attribute_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.attribute_dict[psl_lib.remove_alignment_number(aln_id)].geneType for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class TranscriptType(Attribute):
    """
    Creates a column representing the transcript type
    """
    def run(self):
        self.get_attribute_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.attribute_dict[psl_lib.remove_alignment_number(aln_id)].transcriptType for aln_id
                        in self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class SourceChrom(Attribute):
    """
    Creates a column representing the source chromosome
    """
    def run(self):
        self.get_annotation_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.annotation_dict[psl_lib.remove_alignment_number(aln_id)].chromosome for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class SourceStart(Attribute):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    def run(self):
        self.get_annotation_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.annotation_dict[psl_lib.remove_alignment_number(aln_id)].start for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class SourceStop(Attribute):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    def run(self):
        self.get_annotation_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: self.annotation_dict[psl_lib.remove_alignment_number(aln_id)].stop for aln_id in
                        self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class SourceStrand(Attribute):
    """
    Creates a column representing the source genomic strand.
    """
    def run(self):
        self.get_annotation_dict()
        self.get_alignment_dict()
        results_dict = {aln_id: seq_lib.convert_strand(
                        self.annotation_dict[psl_lib.remove_alignment_number(aln_id)].strand)
                        for aln_id in self.alignment_dict}
        self.dump_attribute_results_to_disk(results_dict)


class DestChrom(Attribute):
    """
    Creates a column representing the dest chromosome
    """
    def run(self):
        self.get_transcript_dict()
        results_dict = {aln_id: self.transcript_dict[aln_id].chromosome for aln_id in self.transcript_dict}
        self.dump_attribute_results_to_disk(results_dict)


class DestStart(Attribute):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    def run(self):
        self.get_transcript_dict()
        results_dict = {aln_id: self.transcript_dict[aln_id].start for aln_id in self.transcript_dict}
        self.dump_attribute_results_to_disk(results_dict)


class DestStop(Attribute):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    def run(self):
        self.get_transcript_dict()
        results_dict = {aln_id: self.transcript_dict[aln_id].stop for aln_id in self.transcript_dict}
        self.dump_attribute_results_to_disk(results_dict)


class DestStrand(Attribute):
    """
    Creates a column representing the dest genomic strand.
    """
    def run(self):
        self.get_transcript_dict()
        results_dict = {aln_id: seq_lib.convert_strand(self.transcript_dict[aln_id].strand) for aln_id in
                        self.transcript_dict}
        self.dump_attribute_results_to_disk(results_dict)


class AlignmentCoverage(Attribute):
    """
    Calculates alignment coverage:

    (matches + mismatches + repeat matches) / qSize

    Reports the value as a REAL between 0 and 1
    """
    def run(self):
        self.get_alignment_dict()
        results_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            results_dict[aln_id] = format_ratio(aln.matches + aln.mismatches + aln.repmatches, aln.q_size)
        self.dump_attribute_results_to_disk(results_dict)


class AlignmentIdentity(Attribute):
    """
    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1
    """
    def run(self):
        self.get_alignment_dict()
        results_dict = {}
        for aln_id, aln in self.alignment_dict.iteritems():
            results_dict[aln_id] = format_ratio(aln.matches + aln.repmatches,
                                                aln.matches + aln.repmatches + aln.mismatches + aln.q_num_insert)
        self.dump_attribute_results_to_disk(results_dict)
