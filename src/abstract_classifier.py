"""
Base classifier classes used by all of the classifiers.
"""
import os
import cPickle as pickle
import itertools

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

    def __init__(self, genome, aln_psl, fasta, ref_fasta, annotation_gp, gencode_attributes, target_gp, ref_genome,
                 tmp_dir):
        # sanity check
        assert all([genome in x for x in [aln_psl, fasta, target_gp]])
        # initialize the Target
        Target.__init__(self)
        self.genome = genome
        self.ref_genome = ref_genome
        self.alnPsl = aln_psl
        self.fasta_path = fasta
        self.ref_fasta_path = ref_fasta
        self.gencode_attributes = gencode_attributes
        self.target_gp = target_gp
        self.annotation_gp = annotation_gp
        self.tmp_dir = tmp_dir
        # these variables will be initialized by individual classifiers as needed
        self.transcripts = None
        self.transcript_dict = None
        self.ref_fasta = None
        self.fasta = None
        self.psls = None
        self.alignment_dict = None
        self.annotations = None
        self.annotation_dict = None
        self.attribute_dict = None

    def get_transcript_dict(self):
        self.transcripts = seq_lib.get_gene_pred_transcripts(self.target_gp)
        self.transcript_dict = seq_lib.transcript_list_to_dict(self.transcripts)

    def get_ref_fasta(self):
        self.ref_fasta = seq_lib.get_sequence_dict(self.ref_fasta_path)

    def get_fasta(self):
        self.fasta = seq_lib.get_sequence_dict(self.fasta_path)

    def get_alignment_dict(self):
        self.psls = psl_lib.read_psl(self.alnPsl)
        self.alignment_dict = psl_lib.get_psl_dict(self.psls)

    def get_annotation_dict(self):
        self.annotations = seq_lib.get_gene_pred_transcripts(self.annotation_gp)
        self.annotation_dict = seq_lib.transcript_list_to_dict(self.annotations)

    @property
    def column(self):
        return self.__class__.__name__

    def dump_results_to_disk(self, classify_dict, details_dict):
        """
        Dumps a pair of classify/details dicts to disk in the globalTempDir for later merging.
        """
        details_dict = sql_lib.collapse_details_dict(details_dict)
        for db, this_dict in itertools.izip(*[["details", "classify"], [details_dict, classify_dict]]):
            base_p = os.path.join(self.tmp_dir, db)
            mkdir_p(base_p)
            p = os.path.join(base_p, self.column)
            with open(p, "wb") as outf:
                pickle.dump(this_dict, outf)


class AbstractAugustusClassifier(AbstractClassifier):
    """
    Overwrites AbstractClassifier to add the extra genePred information and a way to load it.
    """
    def __init__(self, genome, aln_psl, fasta, ref_fasta, annotation_gp, gencode_attributes, target_gp, ref_genome,
                 tmp_dir, augustus_gp):
        AbstractClassifier.__init__(self, genome, aln_psl, fasta, ref_fasta, annotation_gp, gencode_attributes,
                                    target_gp, ref_genome)
        assert self.genome in augustus_gp
        self.augustus_gp = augustus_gp
        self.augustus_transcripts = None
        self.augustus_transcript_dict = None

    def get_augustus_transcript_dict(self):
        self.augustus_transcripts = seq_lib.get_gene_pred_transcripts(self.augustus_gp)
        self.augustus_transcript_dict = seq_lib.transcript_list_to_dict(self.augustus_transcripts)


class Attribute(AbstractClassifier):
    """Need to overwrite the dumpValueDict method for attributes"""
    def get_attribute_dict(self):
        self.attribute_dict = seq_lib.get_transcript_attribute_dict(self.gencode_attributes)

    def dump_attribute_results_to_disk(self, attribute_dict):
        """
        Dumps a attribute dict.
        """
        db = "attributes"
        base_p = os.path.join(self.tmp_dir, db)
        mkdir_p(base_p)
        p = os.path.join(base_p, self.column)
        with open(p, "wb") as outf:
            pickle.dump(attribute_dict, outf)