import re
from itertools import izip
from collections import defaultdict

from lib.general_lib import formatRatio
from src.abstractClassifier import AbstractClassifier

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib


class CodingInsertions(AbstractClassifier):
    """

    Does the alignment introduce insertions that are not a multiple of 3 to the target genome?

    classify mode: reports 1 (TRUE) if so, 0 (FALSE) otherwise

    detail mode: reports BED-format record of each location this happens in.

    Target insertion:

    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def analyzeExons(self, transcript, aln, mult3=False):
        """
        Analyze a Transcript object for coding insertions.
        if mult3 is True, only multiple of 3 insertions are reported.
        """
        insertFlag = False
        tmp = []
        for exon in transcript.exons:
            for i in xrange(exon.start, exon.stop):
                #we have found an insertion
                if insertFlag is False and aln.queryCoordinateToTarget(i) == None:
                    insertSize = 1
                    insertFlag = True
                #insertion continues
                elif insertFlag is True and aln.queryCoordinateToTarget(i) == None:
                    insertSize += 1
                #exiting insertion
                elif insertFlag is True and aln.queryCoordinateToTarget(i) == None:
                    if (annotatedTranscript.transcriptCoordinateToCds(i) is not None or 
                            annotatedTranscript.transcriptCoordinateToCds(i + insertSize) is not None):
                        if insertSize % 3 == 0 and mult3 == True:
                            tmp.append([i - 1, i + insertSize - 1])
                        elif insertSize % 3 != 0 and mult3 == False:
                            tmp.append([i - 1, i + insertSize - 1])
                    insertSize = 0
                    insertFlag = False
        if len(tmp) == 0:
            return None
        else:
            return tmp

    def classifyEntryIter(self, valueDict):
        for aId, entry in valueDict.iteritems():
            if entry is None:
                yield [aId, "0"]
            else:
                yield [aId, "1"]

    def bedEntryIter(self, valueDict):
        tmp = []
        for aId, positions in valueDict.iteritems():
            t = self.transcriptDict[aId]
            chrom = t.chromosomeInterval.chromosome
            for chromStart, chromStop in positions:
                tmp.append([chrom, chromStart, chromStop, psl_lib.removeAlignmentNumber(aId),
                        0, seq_lib.convertStrand(t.strand)])
        return self.detailsEntryIter(tmp)

    def run(self, mult3=False):
        self.getAnnotationDict()
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = defaultdict(list)
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            #annotated transcript coordinates are the same as query coordinates (they are the query)
            annotatedTranscript = self.annotationDict[psl_lib.removeAlignmentNumber(aId)]
            valueDict[aId] = self.analyzeExons(annotatedTranscript, aln, mult3)

        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            if self.details is False:
                sql_lib.updateRows(cur, self.genome, self.getColumn(), self.classifyEntryIter(valueDict))
            else:
                sql_lib.updateRows(cur, self.genome, self.getColumn(), self.bedEntryIter(valueDict))


class CodingMult3Insertions(CodingInsertions):
    """

    See CodingInsertions. Reports all cases where there are multiple of 3 insertions.

    """

    def run(self):
        CodingInsertions.run(self, mult3=True)


class CodingDeletions(AbstractClassifier):
    """

    Does the alignment introduce deletions that are not a multiple of 3 to the target genome?

    classify mode: reports 1 (TRUE) if so, 0 (FALSE) otherwise

    detail mode: reports BED-format record of each location this happens in.
        This BED record is the coordinates of the start-stop, which are always 1 apart.

    Example: (BED entry would say <chrom> 5 6 <name>)

    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA
             012345  67
    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def analyzeExons(self, t, aln, mult3=False):
        tmp = []
        delFlag = False
        for exon in t.exons:
            for i in xrange(exon.start, exon.stop):
                chrom_i = t.transcriptCoordinateToChromosome(i)
                #entering deletion
                if delFlag is False and aln.targetCoordinateToQuery(chrom_i) is None:
                    delSize = 1
                    delFlag = True
                #continuing deletion
                elif delFlag is True and aln.targetCoordinateToQuery(chrom_i) is None:
                    delSize += 1
                #exiting deletion
                elif delFlag is True and aln.targetCoordinateToQuery(chrom_i) is not None:
                    if t.chromosomeCoordinateToCds(chrom_i) is not None:
                        if delSize % 3 == 0 and mult3 is True:
                            tmp.append(i)
                        elif delSize % 3 != 0 and mult3 is False:
                            tmp.append(i)
        if len(tmp) == 0:
            return None
        else:
            return tmp

    def classifyEntryIter(self, valueDict):
        for aId, entry in valueDict.iteritems():
            if entry is None:
                yield ["0", aId]
            else:
                yield ["1", aId]

    def bedEntryIter(self, valueDict):
        tmp = []
        for aId, positions in valueDict.iteritems():
            t = self.transcriptDict[aId]
            chrom = t.chromosomeInterval.chromosome
            if positions is not None:
                for i in positions:
                    chromStart = t.transcriptCoordinateToChromosome(i - 1)
                    chromStop = t.transcriptCoordinateToChromosome(i)
                    tmp.append([aId, [chrom, chromStart, chromStop, aId, 0,
                            seq_lib.convertStrand(t.strand)]])
            else:
                tmp.append([aId, None])
        return self.detailsEntryIter(tmp)

    def run(self, mult3=False):
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = defaultdict(list)
        for aId, aln in self.alignmentDict.iteritems():
            if aId not in self.transcriptDict:
                continue
            transcript = self.transcriptDict[aId]
            valueDict[aId] = self.analyzeExons(transcript, aln, mult3)

        with sql_lib.ExclusiveSqlConnection(self.db) as cur:
            if self.details is False:
                sql_lib.updateRows(cur, self.genome, self.getColumn(), self.classifyEntryIter(valueDict))
            else:
                sql_lib.updateRows(cur, self.genome, self.getColumn(), self.bedEntryIter(valueDict))


class CodingMult3Deletions(CodingDeletions):
    """

    See CodingDeletions. Reports all cases where there are multiple of 3 insertions.

    """

    def run(self):
        CodingDeletions.run(self, mult3=True)


class AlignmentAbutsLeft(AbstractClassifier):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....

    This classifier does not have a details mode.
    Entries are either 1 (TRUE) or 0 (FALSE)

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
                valueDict[aId] = 1
            elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class AlignmentAbutsRight(AbstractClassifier):
    """

    Does the alignment extend off the 3' end of a scaffold?
    (regardless of transcript orientation)

    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|

    This classifier does not have a details mode.
    Entries are either 1 (TRUE) or 0 (FALSE)

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aln.strand == "+" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
                valueDict[aId] = 1
            elif aln.strand == "-" and aln.tStart == 0 and aln.qStart != 0:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class AlignmentCoverage(AbstractClassifier):
    """

    Calculates alignment coverage:

    (matches + mismatches) / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    This classifier does not have a details mode.

    """

    @staticmethod
    def _getClassifierType():
        return "REAL"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches + aln.misMatches, aln.matches + aln.misMatches 
                    + aln.qNumInsert)
        self.simpleUpdateWrapper(valueDict)


class AlignmentIdentity(AbstractClassifier):
    """

    Calculates alignment identity:

    matches / (matches + mismatches + query_insertions)

    Reports the value as a REAL between 0 and 1

    This classifier does not have a details mode.

    """

    @staticmethod
    def _getClassifierType():
        return "REAL"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            valueDict[aId] = formatRatio(aln.matches, aln.matches + aln.misMatches + aln.qNumInsert)
        self.simpleUpdateWrapper(valueDict)


class AlignmentPartialMap(AbstractClassifier):
    """

    Does the query sequence NOT map entirely?

    a.qSize != a.qEnd - a.qStart

    Reports 1 if TRUE and 0 if FALSE

    This classifier does not have a details mode.

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAlignmentDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():
            if aln.qSize != aln.qEnd - aln.qStart:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class BadFrame(AbstractClassifier):
    """

    Looks for CDS sequences that are not a multiple of 3

    classify mode: Reports 1 if TRUE and 0 if FALSE

    details mode: Reports a BED-format record of each CDS which has a bad frame

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self):
        self.getAlignmentDict()
        self.getTranscriptDict()
        valueDict = {}
        for aId, aln in self.alignmentDict.iteritems():        
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.getCdsLength() % 3 != 0:
                if self.details is True:
                    valueDict[aId] = t.getBed()
                else:
                    valueDict[aId] = 1
            else:
                if self.details is False:
                    valueDict[aId] = 0
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)


class TranscriptID(AbstractClassifier):
    """
    Creates a column representing the transcript ID
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        valueDict = {aId : psl_lib.removeAlignmentNumber(aId) for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)
        

class GeneID(AbstractClassifier):
    """
    Creates a column representing the gene ID
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        valueDict = {aId : self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneId 
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class GeneName(AbstractClassifier):
    """
    Creates a column representing the gene name
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        valueDict = {aId : self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneName 
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class GeneType(AbstractClassifier):
    """
    Creates a column representing the gene type
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        valueDict = {aId : self.attributeDict[psl_lib.removeAlignmentNumber(aId)].geneType 
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class TranscriptType(AbstractClassifier):
    """
    Creates a column representing the transcript type
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getAttributeDict()
        valueDict = {aId : self.attributeDict[psl_lib.removeAlignmentNumber(aId)].transcriptType 
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class BeginStart(AbstractClassifier):
    """

    Does the annotated CDS have a start codon (ATG) in the first 3 bases?

    classify mode: Returns 1 if TRUE 0 if FALSE

    details mode: Returns a BED record of the first 3 bases if FALSE

    Value will be NULL if there is unsufficient information, which is defined as:
        1) thickStart == thickStop == 0 (no CDS)
        2) thickStop - thickStart < 3: (no useful CDS annotation)
        3) this alignment was not trans-mapped

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self):
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.thickStart == t.thickStop == 0 or t.thickStop - t.thickStart < 3:
                continue
            s = t.getCds(self.seqDict)
            if s.startswith("ATG"):
                if self.details is False:
                    valueDict[aId] = 1
            else:
                if self.details is False:
                    valueDict[aId] = 0
                else:
                    valueDict[aId].append([t.chromosomeInterval.chromosome, 
                            t.transcriptCoordinateToChromosome(0), 
                            t.transcriptCoordinateToChromosome(3), 
                            aId, 0, seq_lib.convertStrand(t.chromosomeInterval.strand)])
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)


class CdsGap(AbstractClassifier):
    """

    Are any of the CDS introns too short? Too short default is 30 bases.

    classify mode: Returns 1 if TRUE 0 if FALSE

    details mode: Returns BED records of the target genome where these short introns are.

    If mult3 is true, will only report on multiple of 3 gaps.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def detailMode(self, t, mult3, shortIntronSize):
        tmp = []

        for i in xrange(len(t.intronIntervals)):
            #is this intron coding?
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) <= shortIntronSize:
                    if mult3 is True and len(t.intronIntervals[i]) % 3 == 0:
                        for interval in t.intronIntervals:
                            tmp.append(self.parseInterval(interval, t))
                    elif mult3 is False:
                        for interval in t.intronIntervals:
                            tmp.append(self.parseInterval(interval, t))
        if len(tmp) == 0:
            return None
        else:
            return tmp

    def classifyMode(self, t, mult3, shortIntronSize):
        for i in xrange(len(t.intronIntervals)):
            #is this intron coding?
            if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                if len(t.intronIntervals[i]) <= shortIntronSize:
                    if mult3 is True and len(t.intronIntervals[i]) % 3 == 0:
                        return 1
                    elif mult3 is False:
                        return 1
        return 0

    def run(self, mult3=False, shortIntronSize=30):
        self.getTranscriptDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            elif self.details is True:
                valueDict[aId] = self.detailMode(self.transcriptDict[aId], mult3, shortIntronSize)
            else:
                valueDict[aId] = self.classifyMode(self.transcriptDict[aId], mult3, shortIntronSize)
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)
        

class CdsMult3Gap(CdsGap):
    """

    See CdsGap for details. Runs it in mult3 mode.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self, mult3=True, shortIntronSize=30):
        CdsGap.run(self, mult3, shortIntronSize)


class CdsNonCanonSplice(AbstractClassifier):
    """

    Are any of the CDS introns splice sites not of the canonical form
    GT..AG

    classify mode: reports 1 if TRUE, 0 if FALSE

    details mode: reports BED records of each intron which is not canonical

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def badSplice(self, donor, acceptor):
        m = {"GT":"AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue          
            t = self.transcriptDict[aId]
            for i in xrange(len(t.intronIntervals)):
                if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                    if len(t.intronIntervals[i]) > shortIntronSize:
                        donor = self.seqDict[t.chromosomeInterval.chromosome][t.intronIntervals[i].start : t.intronIntervals[i].start + 2]    
                        acceptor = self.seqDict[t.chromosomeInterval.chromosome][t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                        if self.details is True and self.badSplice(donor, acceptor) is True:
                            valueDict[aId].append(self.parseInterval(t.intronIntervals[i], t))
                        elif self.details is False and self.badSplice(donor, acceptor) is True:
                            valueDict[aId] = 1
                            break
            if self.details is False and aId not in valueDict:
                valueDict[aId] = 0
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)    


class CdsUnknownSplice(CdsNonCanonSplice):
    """

    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    subclasses cdsNonCanonSplice and just replaces the badSplice function

    classify mode: reports 1 if TRUE, 0 if FALSE

    details mode: reports BED records of each intron which is not canonical

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def badSplice(self, donor, acceptor):
        m = {"GT":"AG", "GC":"AG", "AT":"AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        CdsNonCanonSplice.run(self)


class EndStop(AbstractClassifier):
    """

    Looks at the end of the coding region (thickEnd) and sees if the last
    three bases are a stop codon ('TAA', 'TGA', 'TAG')

    classify mode: Returns 1 if TRUE 0 if FALSE

    details mode: Reports a BED record for the last 3 bases if FALSE

    Value will be NULL if there is unsufficient information, which is defined as:
        1) thickStop - thickStart < 3: (no useful CDS annotation)
        2) this alignment was not trans-mapped

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self):
        stopCodons = ('TAA', 'TGA', 'TAG')
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.thickStop - t.thickStart < 3:
                continue
            s = t.getCds(self.seqDict)[-3:]
            if s in stopCodons:
                if self.details is False:
                    valueDict[aId] = 1
            else:
                if self.details is False:
                    valueDict[aId] = 0
                else:
                    l = t.getCdsLength()
                    valueDict[aId] = [t.chromosomeInterval.chromosome, 
                            t.cdsCoordinateToChromosome(l-3), t.cdsCoordinateToChromosome(l), 
                            aId, 0, seq_lib.convertStrand(t.chromosomeInterval.strand)]
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)

class InFrameStop(AbstractClassifier):
    """

    Reports on in frame stop codons for each transcript.

    classify mode: Reports 1 if TRUE (has in frame stop), 0 if FALSE

    details mode: Reports a BED record for the first in frame stop

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self):
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            #make sure this transcript has CDS
            #and more than 2 codons - can't have in frame stop without that
            cds_size = t.getCdsLength()
            if cds_size >= 9:
                for i in xrange(9, cds_size - 3, 3):
                    c = t.cdsCoordinateToAminoAcid(i, self.seqDict)
                    if c == "*":
                        if self.details is True:
                            valueDict[aId] = [t.chromosomeInterval.chromosome, 
                                    t.cdsCoordinateToChromosome(i),
                                    t.cdsCoordinateToChromosome(i+3), 
                                    aId, 0, seq_lib.convertStrand(t.chromosomeInterval.strand)]
                        else:
                            valueDict[aId] = 1
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)


class NoCds(AbstractClassifier):
    """

    Looks to see if this transcript actually has a CDS, which is defined as having a
    thickStop-thickStart region of at least 1 codon. Adjusting cdsCutoff can change this.

    This classifier has no details mode.

    Reports a 1 if TRUE, 0 if FALSE.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self, cdsCutoff=3):
        self.getTranscriptDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if t.getCdsLength() < cdsCutoff:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class MinimumCdsSize(NoCds):
    """

    The smallest ORFs in any species are >10AA. So, we will flag any CDS smaller than this.

    Inherits NoCds and modifies cdsCutoff to do this.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        NoCds.run(self, cdsCutoff=30)


class ScaffoldGap(AbstractClassifier):
    """

    Does this alignment span a scaffold gap? (Defined as a 100bp run of Ns)

    This classifier does not have a details mode.

    Reports 1 if TRUE, 0 if FALSE

    """

    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAlignmentDict()
        self.getSeqDict()
        valueDict = {}
        r = re.compile("[N]{100}")
        for aId, aln in self.alignmentDict.iteritems():
            destSeq = self.seqDict[aln.tName][aln.tStart : aln.tEnd].upper()
            if re.search(r, destSeq) is not None:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class SourceChrom(AbstractClassifier):
    """
    Creates a column representing the source chromosome
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        valueDict = {aId : self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.chromosome
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)

class SourceStart(AbstractClassifier):
    """
    Creates a column representing the source genomic start location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        valueDict = {aId : self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.start
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class SourceStop(AbstractClassifier):
    """
    Creates a column representing the source genomic stop location.
    (+) strand value, so always smaller than sourceEnd.
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getAnnotationDict()
        valueDict = {aId : self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.stop
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class SourceStrand(AbstractClassifier):
    """
    Creates a column representing the source genomic strand.
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getAnnotationDict()
        valueDict = {aId : seq_lib.convertStrand(self.annotationDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.strand)
                for aId in self.aIds}
        self.simpleUpdateWrapper(valueDict)


class DestChrom(AbstractClassifier):
    """
    Creates a column representing the dest chromosome
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId : self.transcriptDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.chrom
                for aId in self.aIds if psl_lib.removeAlignmentNumber(aId) in self.transcriptDict}
        self.simpleUpdateWrapper(valueDict)


class DestStart(AbstractClassifier):
    """
    Creates a column representing the dest genomic start location.
    (+) strand value, so always smaller than destEnd.
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId : self.transcriptDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.start
                for aId in self.aIds if psl_lib.removeAlignmentNumber(aId) in self.transcriptDict}
        self.simpleUpdateWrapper(valueDict)


class DestStop(AbstractClassifier):
    """
    Creates a column representing the dest genomic stop location.
    (+) strand value, so always larger tha destStart
    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId : self.transcriptDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.stop
                for aId in self.aIds if psl_lib.removeAlignmentNumber(aId) in self.transcriptDict}
        self.simpleUpdateWrapper(valueDict)


class DestStrand(AbstractClassifier):
    """
    Creates a column representing the dest genomic strand.
    """
    @staticmethod
    def _getClassifierType():
        return "TEXT"

    def run(self):
        self.getTranscriptDict()
        valueDict = {aId : seq_lib.convertStrand(self.transcriptDict[psl_lib.removeAlignmentNumber(aId)].chromosomeInterval.strand)
                for aId in self.aIds if psl_lib.removeAlignmentNumber(aId) in self.transcriptDict}
        self.simpleUpdateWrapper(valueDict)


class UnknownBases(AbstractClassifier):
    """

    Does this alignment contain Ns in the target genome?

    Classify mode: Reports 1 if TRUE, 0 if FALSE

    This classifier does not have a details mode.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self, cds=False):
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = {}
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            if cds is True:
                s = t.getCds(self.seqDict)
            else:
                s = t.getMRna(self.seqDict)
            if "N" in s:
                valueDict[aId] = 1
            else:
                valueDict[aId] = 0
        self.simpleUpdateWrapper(valueDict)


class UnknownCdsBases(UnknownBases):
    """

    Inherits Unknown Bases and sets the cds flag to True.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    def run(self):
        UnknownBases.run(self, cds=True)


class UtrGap(AbstractClassifier):
    """

    Are any UTR introns too short?

    Classify mode: Reports 1 if TRUE, 0 if FALSE

    Details mode: Reports BED records of each UTR intron that is too short.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def run(self, shortIntronSize=30):
        self.getTranscriptDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue
            t = self.transcriptDict[aId]
            for i in xrange(len(t.intronIntervals)):
                if t.exons[i].containsCds() is False and t.exons[i+1].containsCds() is False:
                    if len(t.intronIntervals[i]) <= shortIntronSize:
                        if self.details is False:
                            valueDict[aId] = 1
                            break
                        else:
                            valueDict[aId].append(self.parseInterval(t.intronIntervals[i], t))
            if aId not in valueDict:
                if self.details is False:
                    valueDict[aId] = 0
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)


class UtrNonCanonSplice(AbstractClassifier):
    """

    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    classify mode: reports 1 if TRUE, 0 if FALSE

    details mode: reports BED records of each intron which is not canonical

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    TODO: this class is identical to the CDS version, replacing two True with False.

    Probably should merge them.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def badSplice(self, donor, acceptor):
        m = {"GT":"AG"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        self.getTranscriptDict()
        self.getSeqDict()
        valueDict = defaultdict(list)
        for aId in self.aIds:
            if aId not in self.transcriptDict:
                continue          
            t = self.transcriptDict[aId]
            for i in xrange(len(t.intronIntervals)):
                if t.exons[i].containsCds() is False and t.exons[i+1].containsCds() is False:
                    if len(t.intronIntervals[i]) > shortIntronSize:
                        donor = self.seqDict[t.chromosomeInterval.chromosome][t.intronIntervals[i].start : t.intronIntervals[i].start + 2]    
                        acceptor = self.seqDict[t.chromosomeInterval.chromosome][t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                        if self.details is True and self.badSplice(donor, acceptor) is True:
                            valueDict[aId].append(self.parseInterval(t.intronIntervals[i], t))
                        elif self.details is False and self.badSplice(donor, acceptor) is True:
                            valueDict[aId] = 1
                            break
            if self.details is False and aId not in valueDict:
                valueDict[aId] = 0
        if self.details is False:
            self.simpleUpdateWrapper(valueDict)
        else:
            self.simpleBedUpdateWrapper(valueDict)  


class UtrUnknownSplice(UtrNonCanonSplice):
    """

    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    subclasses cdsNonCanonSplice and just replaces the badSplice function

    classify mode: reports 1 if TRUE, 0 if FALSE

    details mode: reports BED records of each intron which is not canonical

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    """
    @staticmethod
    def _getClassifierType():
        return "INTEGER"

    @staticmethod
    def _getDetailsType():
        return "TEXT"

    def badSplice(self, donor, acceptor):
        m = {"GT":"AG", "GC":"AG", "AT":"AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, shortIntronSize=30):
        UtrNonCanonSplice.run(self)