import os
import cPickle as pickle

from jobTree.scriptTree.target import Target

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class AbstractClassifier(Target):
    def __init__(self, genome, alnPsl, fasta, refSeqTwoBit, annotationBed, gencodeAttributeMap, geneCheckBed, refGenome,
                 primaryKey, outDir):
        #initialize the Target
        Target.__init__(self)
        #primary key this will be keyed on (alignmentID usually)
        self.primaryKey = primaryKey
        self.genome = genome
        self.refGenome = refGenome
        self.alnPsl = alnPsl
        self.fasta = fasta
        self.refSeqTwoBit = refSeqTwoBit
        self.gencodeAttributeMap = gencodeAttributeMap
        self.geneCheckBed = geneCheckBed
        self.annotationBed = annotationBed
        self.outDir = os.path.join(outDir, self.genome)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

        # alignment IDs
        self.aIds = set(x.split()[9] for x in open(alnPsl))

        # for details classifiers, color codes for types of names for BED record
        self.colors = {'input': '219,220,222',     # grey
                       'mutation': '132,35,27',    # red-ish
                       'assembly': '167,206,226',  # light blue
                       'alignment': '35,125,191',  # blue
                       'synon': '163,116,87',      # light brown
                       'nonsynon': '181,216,139',  # avocado
                       'generic': '152,156,45'     # grey-yellow
                      }

    def getTranscriptDict(self):
        self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)
        self.transcriptDict = seq_lib.transcriptListToDict(self.transcripts, noDuplicates=True)

    def getRefDict(self):
        self.refDict = seq_lib.readTwoBit(self.refSeqTwoBit)

    def getSeqDict(self):
        self.seqDict = seq_lib.getSequenceDict(self.fasta)

    def getAlignmentDict(self):
        self.psls = psl_lib.readPsl(self.alnPsl)
        self.alignmentDict = psl_lib.getPslDict(self.psls, noDuplicates=True)

    def getAnnotationDict(self):
        self.annotations = seq_lib.getTranscripts(self.annotationBed)
        self.annotationDict = seq_lib.transcriptListToDict(self.annotations, noDuplicates=True)

    def getColumn(self):
        return self.__class__.__name__

    def dumpValueDicts(self, classifyDict, detailsDict):
        """
        Dumps a pair of classify/details dicts to disk in the globalTempDir for later merging.
        """
        #with open(os.path.join(self.getGlobalTempDir(), "Details" + self.getColumn() + self.genome), "wb") as outf:
        with open(os.path.join(self.outDir, "Details" + self.getColumn() + self.genome), "wb") as outf:
            pickle.dump(detailsDict, outf)
        #with open(os.path.join(self.getGlobalTempDir(), "Classify" + self.getColumn() + self.genome), "wb") as outf:
        with open(os.path.join(self.outDir, "Classify" + self.getColumn() + self.genome), "wb") as outf:
            pickle.dump(classifyDict, outf)            


class Attribute(AbstractClassifier):
    """Need to overwrite the dumpValueDict method for attributes"""
    def getAttributeDict(self):
        self.attributeDict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)        

    def dumpValueDict(self, valueDict):
        """
        Dumps a attribute dict.
        """
        #with open(os.path.join(self.getGlobalTempDir(), "Attribute" + self.genome), "wb") as outf:
        with open(os.path.join(self.outDir, "Attribute" + self.getColumn() + self.genome), "wb") as outf:
            pickle.dump(valueDict, outf)
