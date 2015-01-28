from jobTree.scriptTree.target import Target

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class AbstractClassifier(Target):
    def __init__(self, genome, alnPsl, seqTwoBit, refSeqTwoBit, annotationBed,   
                gencodeAttributeMap, geneCheckBed, outDb, refGenome, primaryKey, details=False):
        #initialize the Target
        Target.__init__(self)

        #flag to turn on details mode for classifiers that have it
        self.details = details

        #primary key this will be keyed on (alignmentID usually)
        self.primaryKey = primaryKey

        self.genome = genome
        self.refGenome = refGenome
        self.alnPsl = alnPsl
        self.seqTwoBit = seqTwoBit
        self.refSeqTwoBit = refSeqTwoBit
        self.gencodeAttributeMap = gencodeAttributeMap
        self.geneCheckBed = geneCheckBed
        self.annotationBed = annotationBed
        self.db = outDb

        #alignment IDs
        self.aIds = alnIds = set(x.split()[9] for x in open(alnPsl))


    def getAttributeDict(self):
        self.attributeDict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)

    def getTranscriptDict(self):
        self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)
        self.transcriptDict = seq_lib.transcriptListToDict(self.transcripts, noDuplicates=True)

    def getSeqDict(self):
        self.seqDict = seq_lib.readTwoBit(self.seqTwoBit)

    def getAlignmentDict(self):
        self.psls = psl_lib.readPsl(self.alnPsl)
        self.alignmentDict = psl_lib.getPslDict(self.psls, noDuplicates=True)

    def getAnnotationDict(self):
        self.annotations = seq_lib.getTranscripts(self.annotationBed)
        self.annotationDict = seq_lib.transcriptListToDict(self.annotations, noDuplicates=True)

    def detailsEntryIter(self, valueIter):
        """
        General case for converting a valueIter to a string entry representing 1 or more BED records
        valueIters are generally lists of lists where each sublist represents a BED record
        """
        for aId, entry in valueIter:
            if entry is None:
                yield [None, aId]
            elif type(entry[0]) != list:
                #only one entry
                bedEntries = ["\t".join(map(str,entry))]
            else:
                bedEntries = ["\t".join(map(str, x)) for x in entry]
                yield ["\n".join(bedEntries), aId]

    def getColumn(self):
        return self.__class__.__name__

    def invertDict(self, d):
        for x in d.iteritems():
            yield x[::-1]

    def simpleUpdateWrapper(self, valueDict):
        """
        If your classifier is going to do a simple 1-1 update with a valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            sql_lib.updateRows(cur, self.genome, self.primaryKey, self.getColumn(), 
                    self.invertDict(valueDict))

    def simpleBedUpdateWrapper(self, valueDict):
        """
        If your details-mode classifier has BED records for values in its valueDict, use this.
        """
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            sql_lib.updateRows(cur, self.genome, self.primaryKey, self.getColumn(), 
                    self.detailsEntryIter(valueDict.iteritems()))     

    def parseInterval(self, interval, t):
        """
        If you are turning interval objects into BED records, look here. t is a transcript object.
        """
        return [interval.chromosome, interval.start, interval.stop, 
                t.name, 0, seq_lib.convertStrand(interval.strand)]

    def transcriptToBed(self, t):
        """
        Convenience function for pulling BED tokens from a Transcript object.
        """
        return t.getBed()