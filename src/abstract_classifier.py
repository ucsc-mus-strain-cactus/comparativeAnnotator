import os

from jobTree.scriptTree.target import Target
from sonLib.bioio import logger
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class AbstractClassifier(Target):
    def __init__(self, genome, alnPsl, seqFasta, annotationBed, gencodeAttributeMap,  
                geneCheckBed, outDir, refGenome, primaryKey):
        #initialize the Target
        Target.__init__(self)

        #store basic information
        self.genome = genome
        self.refGenome = refGenome
        self.alnPsl = alnPsl
        self.seqFasta = seqFasta
        self.annotationBed = annotationBed
        self.gencodeAttributeMap = gencodeAttributeMap
        self.geneCheckBed = geneCheckBed
        self.primary_key = primaryKey
        self.db = os.path.join(outDir, self.genome + ".db")

    def get_alignment_ids(self):
        self.alignment_ids = set(x.split()[9] for x in open(self.alnPsl))

    def get_original_transcripts(self):
        self.original_transcripts = seq_lib.getTranscripts(self.annotationBed)

    def get_transcript_attributes(self):
        self.attribute_dict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)

    def get_original_transcript_dict(self):
        if not hasattr(self, 'original_transcripts'):
            self.original_transcripts = seq_lib.getTranscripts(self.annotationBed)
        self.original_transcript_dict = seq_lib.transcriptListToDict(self.original_transcripts, noDuplicates=True)

    def get_transcripts(self):
        self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)

    def get_transcript_dict(self):
        if not hasattr(self, 'transcripts'):
            self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)
        self.transcript_dict = seq_lib.transcriptListToDict(self.transcripts, noDuplicates=True)

    def get_seq_dict(self):
        self.seq_dict = seq_lib.readTwoBit(self.seqFasta)

    def get_alignments(self):
        self.alignments = psl_lib.readPsl(self.alnPsl)

    def get_alignment_dict(self):
        if not hasattr(self, 'alignments'):
            self.get_alignments()
        self.alignment_dict = psl_lib.getPslDict(self.alignments, noDuplicates=True)

    def upsert_wrapper(self, alignmentName, value):
        """convenience wrapper for upserting into a column in the sql lib.
        So you don't have to call __name__, self.primaryKey, etc each time"""
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            sql_lib.upsert(cur, self.genome, self.primary_key, alignmentName, 
                    self.__class__.__name__, str(value))

    def upsert_dict_wrapper(self, d):
        """even more convenient wrapper for upserting. Assumes input is a dict 
        mapping alignment names to a value.
        """
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            for aln, value in d.iteritems():
                sql_lib.upsert(cur, self.genome, self.primary_key, aln,
                        self.__class__.__name__, str(value))