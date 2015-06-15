import os
import cPickle as pickle

from jobTree.scriptTree.target import Target

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib


class AbstractClassifier(Target):
    colors = {'input': '219,220,222',     # grey
               'mutation': '132,35,27',    # red-ish
               'assembly': '167,206,226',  # light blue
               'alignment': '35,125,191',  # blue
               'synon': '163,116,87',      # light brown
               'nonsynon': '181,216,139',  # avocado
               'generic': '152,156,45'     # grey-yellow
              }

    def __init__(self, genome, alnPsl, fasta, refFasta, annotationGp, gencodeAttributeMap, targetGp, refGenome,
                 primaryKey, outDir):
        # initialize the Target
        Target.__init__(self, memory=16 * 1024 * 1024 * 1024)  # 16GB RAM per job
        # primary key this will be keyed on (AlignmentId usually)
        self.primaryKey = primaryKey
        self.genome = genome
        self.refGenome = refGenome
        self.alnPsl = alnPsl
        self.fasta = fasta
        self.refFasta = refFasta
        self.gencodeAttributeMap = gencodeAttributeMap
        self.targetGp = targetGp
        self.annotationGp = annotationGp
        self.outDir = os.path.join(outDir, self.genome)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def getTranscriptDict(self):
        self.transcripts = seq_lib.getGenePredTranscripts(self.targetGp)
        self.transcriptDict = seq_lib.transcriptListToDict(self.transcripts, noDuplicates=True)

    def getRefDict(self):
        self.refDict = seq_lib.getSequenceDict(self.refFasta)

    def getSeqDict(self):
        self.seqDict = seq_lib.getSequenceDict(self.fasta)

    def getAlignmentDict(self):
        self.psls = psl_lib.readPsl(self.alnPsl)
        self.alignmentDict = psl_lib.getPslDict(self.psls, noDuplicates=True)

    def getAnnotationDict(self):
        self.annotations = seq_lib.getGenePredTranscripts(self.annotationGp)
        self.annotationDict = seq_lib.transcriptListToDict(self.annotations, noDuplicates=True)

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


class Attribute(AbstractClassifier):
    """Need to overwrite the dumpValueDict method for attributes"""
    def getAttributeDict(self):
        self.attributeDict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)        

    def dumpValueDict(self, valueDict):
        """
        Dumps a attribute dict.
        """
        with open(os.path.join(self.outDir, "Attribute" + self.column + self.genome), "wb") as outf:
            pickle.dump(valueDict, outf)