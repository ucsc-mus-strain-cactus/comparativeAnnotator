"""
Base classifier classes used by all of the classifiers.
"""
import os
import cPickle as pickle
import itertools
import copy_reg
import types
from collections import defaultdict

from jobTree.scriptTree.target import Target

import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sql_lib as sql_lib
from lib.general_lib import mkdir_p

__author__ = "Ian Fiddes"


class AbstractClassifier(Target):
    colors = {'input': '219,220,222',     # grey
              'mutation': '132,35,27',    # red-ish
              'assembly': '167,206,226',  # light blue
              'alignment': '35,125,191',  # blue
              'synon': '163,116,87',      # light brown
              'nonsynon': '181,216,139',  # avocado
              'generic': '152,156,45'     # grey-yellow
              }

    def __init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir):
        # initialize the Target
        Target.__init__(self)
        self.ref_genome = ref_genome
        self.ref_fasta = ref_fasta
        self.annotation_gp = annotation_gp
        self.tmp_dir = tmp_dir
        # these variables will be initialized once the jobs have begun to not pickle all of this stuff needlessly
        self.annotation_dict = None
        self.ref_seq_dict = None
        # these dictionaries will be filled out by each classifier
        self.classify_dict = {}
        self.details_dict = defaultdict(list)

    def get_fasta(self):
        self.ref_seq_dict = seq_lib.get_sequence_dict(self.ref_fasta)

    def get_annotation_dict(self):
        self.annotation_dict = seq_lib.get_transcript_dict(self.annotation_gp)

    def annotation_iterator(self):
        """
        Convenience function for iterating over a dictionary of Transcript objects
        """
        if self.annotation_dict is None:
            self.get_annotation_dict()
        for ens_id, a in self.annotation_dict.iteritems():
            yield ens_id, a

    @property
    def column(self):
        return self.__class__.__name__

    def dump_results_to_disk(self):
        """
        Dumps a pair of classify/details dicts to disk in the globalTempDir for later merging.
        """
        details_dict = sql_lib.collapse_details_dict(self.details_dict)
        for db, this_dict in itertools.izip(*[["details", "classify"], [details_dict, self.classify_dict]]):
            base_p = os.path.join(self.tmp_dir, db)
            mkdir_p(base_p)
            p = os.path.join(base_p, self.column)
            with open(p, "wb") as outf:
                pickle.dump(this_dict, outf)


class AbstractAlignmentClassifier(AbstractClassifier):
    """
    Subclasses AbstractClassifier for alignment classifications
    """
    colors = {'input': '219,220,222',     # grey
              'mutation': '132,35,27',    # red-ish
              'assembly': '167,206,226',  # light blue
              'alignment': '35,125,191',  # blue
              'synon': '163,116,87',      # light brown
              'nonsynon': '181,216,139',  # avocado
              'generic': '152,156,45'     # grey-yellow
              }

    def __init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir, tgt_genome, aln_psl, tgt_fasta, tgt_gp):
        AbstractClassifier.__init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir)
        self.genome = tgt_genome
        self.aln_psl = aln_psl
        self.tgt_gp = tgt_gp
        self.tgt_fasta = tgt_fasta
        self.transcript_dict = None
        self.alignment_dict = None
        self.seq_dict = None

    def get_fasta(self):
        self.seq_dict = seq_lib.get_sequence_dict(self.tgt_fasta)
        self.ref_seq_dict = seq_lib.get_sequence_dict(self.ref_fasta)

    def get_alignment_dict(self):
        if self.transcript_dict is None:
            self.get_transcript_dict()
        self.alignment_dict = psl_lib.get_alignment_dict(self.aln_psl, filter=self.transcript_dict.viewkeys())

    def get_transcript_dict(self):
        self.transcript_dict = seq_lib.get_transcript_dict(self.tgt_gp)

    def transcript_iterator(self):
        """
        Convenience function for iterating over a dictionary of Transcript objects
        """
        if self.transcript_dict is None:
            self.get_transcript_dict()
        for aln_id, t in self.transcript_dict.iteritems():
            yield aln_id, t

    def alignment_iterator(self):
        """
        Convenience function for iterating over a dictionary of Transcript objects
        """
        if self.alignment_dict is None:
            self.get_alignment_dict()
        for aln_id, aln in self.alignment_dict.iteritems():
            yield aln_id, aln

    def alignment_transcript_iterator(self):
        """
        Convenience function for iterating over both Transcript and PslRow objects at once
        """
        for aln_id, aln in self.alignment_iterator():
            t = self.transcript_dict[aln_id]
            yield aln_id, aln, t

    def alignment_transcript_annotation_iterator(self):
        """
        Convenience function for iterating over alignment, ref transcript and tgt transcript
        """
        if self.annotation_dict is None:
            self.get_annotation_dict()
        for aln_id, aln, t in self.alignment_transcript_iterator():
            a = self.annotation_dict[psl_lib.remove_alignment_number(aln_id)]
            yield aln_id, aln, t, a


class AbstractAugustusClassifier(AbstractAlignmentClassifier):
    """
    Subclasses AbstractClassifier for Augustus classifications
    """
    def __init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir, tgt_genome, aln_psl, tgt_fasta, tgt_gp,
                 augustus_gp):
        AbstractAlignmentClassifier.__init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir, tgt_genome, aln_psl,
                                             tgt_fasta, tgt_gp)
        self.augustus_gp = augustus_gp
        self.augustus_transcript_dict = None

    def get_augustus_transcript_dict(self):
        self.augustus_transcript_dict = seq_lib.get_transcript_dict(self.augustus_gp)

    def augustus_transcript_iterator(self):
        if self.augustus_transcript_dict is None:
            self.get_augustus_transcript_dict()
        for aug_id, aug_t in self.augustus_transcript_dict.iteritems():
            yield aug_id, aug_t

    def augustus_transcript_transmap_iterator(self):
        if self.transcript_dict is None:
            self.get_transcript_dict()
        for aug_id, aug_t in self.augustus_transcript_iterator():
            t = self.transcript_dict[psl_lib.remove_augustus_alignment_number(aug_id)]
            yield aug_id, aug_t, t


class Attribute(AbstractAlignmentClassifier):
    """
    Subclasses AbstractClassifier to build the Attributes database
    """
    def __init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir, tgt_genome, aln_psl, tgt_fasta, tgt_gp,
                 gencode_attributes):
        AbstractAlignmentClassifier.__init__(self, ref_fasta, annotation_gp, ref_genome, tmp_dir, tgt_genome, aln_psl,
                                             tgt_fasta, tgt_gp)
        self.gencode_attributes = gencode_attributes
        self.attribute_dict = None

    def get_attribute_dict(self):
        self.attribute_dict = seq_lib.get_transcript_attribute_dict(self.gencode_attributes)

    def attribute_iterator(self):
        if self.attribute_dict is None:
            self.get_attribute_dict()
        for ens_id, data in self.attribute_dict.iteritems():
            yield ens_id, data

    def dump_attribute_results_to_disk(self, results_dict):
        """
        Dumps a attribute dict.
        """
        db = "attributes"
        base_p = os.path.join(self.tmp_dir, db)
        mkdir_p(base_p)
        p = os.path.join(base_p, self.column)
        with open(p, "wb") as outf:
            pickle.dump(results_dict, outf)