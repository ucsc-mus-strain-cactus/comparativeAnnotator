"""
Base classifier classes used by all of the classifiers.
"""
import os
import cPickle as pickle

from jobTree.scriptTree.target import Target

import lib.seq_lib as seq_lib
import lib.psl_lib as psl_lib

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
outDir):
        # sanity check
        assert all([genome in x for x in [aln_psl, fasta, target_gp]])
        # initialize the Target
        Target.__init__(self)
        # primary key this will be keyed on (AlignmentId usually)
        self.primaryKey = primaryKey
        self.genome = genome
        self.refGenome = ref_genome
        self.alnPsl = aln_psl
        self.fasta = fasta
        self.refFasta = ref_fasta
        self.gencodeAttributeMap = gencode_attributes
        self.targetGp = target_gp
        self.annotationGp = annotation_gp
        self.outDir = os.path.join(outDir, self.genome)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def getTranscriptDict(self):
        self.transcripts = seq_lib.get_gene_pred_transcripts(self.targetGp)
        self.transcriptDict = seq_lib.transcript_list_to_dict(self.transcripts, noDuplicates=True)

    def getRefDict(self):
        self.refDict = seq_lib.get_sequence_dict(self.refFasta)

    def getSeqDict(self):
        self.seqDict = seq_lib.get_sequence_dict(self.fasta)

    def getAlignmentDict(self):
        self.psls = psl_lib.read_psl(self.alnPsl)
        self.alignmentDict = psl_lib.get_psl_dict(self.psls, noDuplicates=True)

    def getAnnotationDict(self):
        self.annotations = seq_lib.get_gene_pred_transcripts(self.annotationGp)
        self.annotationDict = seq_lib.transcript_list_to_dict(self.annotations, noDuplicates=True)

    @property
    def column(self):
        return self.__class__.__name__

    def dumpValueDicts(self, classifyDict, detailsDict):
        """
        Dumps a pair of classify/details dicts to disk in the globalTempDir for later merging.
        """
        with open(os.path.join(self.outDir, "Details" + self.column + self.genome), "wb") as outf:
            pickle.dump(detailsDict, outf)
        with open(os.path.join(self.outDir, "Classify" + self.column + self.genome), "wb") as outf:
            pickle.dump(classifyDict, outf)


class AbstractAugustusClassifier(AbstractClassifier):
    """
    Overwrites AbstractClassifier to add the extra genePred information and a way to load it.
    """
    def __init__(self, genome, aln_psl, fasta, ref_fasta, annotation_gp, gencode_attributes, target_gp, ref_genome,
                 primaryKey, outDir, augustusGp):
        AbstractClassifier.__init__(self, genome, aln_psl, fasta, ref_fasta, annotation_gp, gencode_attributes, target_gp,
                                    ref_genome, primaryKey, outDir)
        assert self.genome in augustusGp
        self.augustusGp = augustusGp

    def getAugustusTranscriptDict(self):
        self.augustusTranscripts = seq_lib.get_gene_pred_transcripts(self.augustusGp)
        self.augustusTranscriptDict = seq_lib.transcript_list_to_dict(self.augustusTranscripts, noDuplicates=True)


class Attribute(AbstractClassifier):
    """Need to overwrite the dumpValueDict method for attributes"""
    def getAttributeDict(self):
        self.attributeDict = seq_lib.get_transcript_attribute_dict(self.gencodeAttributeMap)

    def dumpValueDict(self, valueDict):
        """
        Dumps a attribute dict.
        """
        with open(os.path.join(self.outDir, "Attribute" + self.column + self.genome), "wb") as outf:
            pickle.dump(valueDict, outf)