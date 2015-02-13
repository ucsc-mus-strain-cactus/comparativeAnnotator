import os

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from src.abstractClassifier import AbstractClassifier


transcripts = seq_lib.getTranscripts("datafiles/Rattus.gene-check.bed")
transcriptDict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
annotations = seq_lib.getTranscripts("datafiles/wgEncodeGencodeBasicVM2.gene-check.bed")
annotationDict = seq_lib.transcriptListToDict(annotations, noDuplicates=True)
alignments = psl_lib.readPsl("datafiles/Rattus.filtered.psl")
alignmentDict = psl_lib.getPslDict(alignments, noDuplicates=True)
seqDict = seq_lib.readTwoBit("datafiles/Rattus.2bit")


mult3 = False
records = []
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    exonStarts = [x.start for x in a.exons]
    prevTargetPos = None
    for query_i in xrange(len(a)):
        if query_i in exonStarts:
            prevTargetPos = None
        target_i = aln.queryCoordinateToTarget(query_i)
        if target_i is None:
            #deletion
            continue
        if prevTargetPos is not None and abs(target_i - prevTargetPos) != 1:
            #insertion
            start = min(prevTargetPos, target_i) + 1
            stop = max(prevTargetPos, target_i)
            t = transcriptDict[aId]           
            if mult3 is True and start - stop % 3 == 0:
                records.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
            else:
                records.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, self.rgb(), self.getColumn()))
        prevTargetPos = target_i

mult3 = False
records = []
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    t = transcriptDict[aId]
    prevQueryPos = None
    for transcript_i in xrange(len(t)):
        target_i = t.transcriptCoordinateToChromosome(transcript_i)
        query_i = aln.targetCoordinateToQuery(target_i)
        if query_i is None:
            continue
        if prevQueryPos is not None and abs(query_i - prevQueryPos) != 0:
            #make sure this isn't just an intron
            break
            #if seq_lib.ChromosomeInterval(aln.qName, min(prevQueryPos, query_i), max(prevQueryPos, query_i), True) not in a.intronIntervals:
            #    break
        prevQueryPos = query_i


def simplePsl(strand, qSize, qStart, qEnd, tSize, tStart, tEnd,
              blockSizes, qStarts, tStarts, qName='query', tName='target'):
    """ Given a few of the fields, create a PslRow object.
    """
    line = ('%d %d %d %d %d %d %d %d %s %s %d %d %d %s %d %d %d %d %s %s %s'% (1, 0, 0, 0, 0, 0, 0, 0, strand, qName, qSize, qStart, qEnd,tName, tSize, tStart, tEnd, len(blockSizes),','.join([str(b) for b in blockSizes]),','.join([str(b) for b in qStarts]),','.join([str(b) for b in tStarts]),))
    return line

def bedLine(chrom, chromStart, chromEnd, name, score=None, strand=None,
            thickStart=None, thickEnd=None, itemRgb=None, blockCount=None,
            blockSizes=None, blockStarts=None):
    """ Give the fields, create a bed line string
    """
    s = ('%s %d %d %s'
       % (chrom, chromStart, chromEnd, name))
    if score is not None:
        for v in [strand, thickStart, thickEnd, itemRgb,
              blockCount, blockSizes, blockStarts]:
            assert(v is not None)
        s += (' %d %s %d %d %s %d %s %s'
          % (score, strand, thickStart, thickEnd, itemRgb, blockCount,
          blockSizes, blockStarts))
    return s

    ##########
    #            0          11
    # ref        ATGATCCAATGA  query
    # exons       ****  ****
    # non ref    ATGATTAA--GA  target
    #            0          9
    #####

    ##########
    #           0 1 2 3 4 5 6 7 8 9 1011121314151617
    # query     G T A T T G G C T T G G A C
    # target    G T A T T - - C T T G G A C C T A A G

aln = psl_lib.PslRow(simplePsl("+", 14, 0, 14, 290094216, 0, 12, [5,7], [0,7], [0,5], qName="query", tName="chr1"))
t = seq_lib.Transcript(bedLine("chr1", 0, 12, "query", 1000, "+", 0, 12, "128,0,0", 1, 12, 0).split())
a = seq_lib.Transcript(bedLine("chr1", 0, 14, "query", 1000, "+", 0, 14, "128,0,0", 1, 14, 0).split())

mult3 = False
records = []
count = 0
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    t = transcriptDict[aId]        
    records = []
    deleteFlag = False
    for query_i in xrange(len(a)):
        target_i = aln.queryCoordinateToTarget(query_i)
        if target_i is None and deleteFlag is False:
            #entering deletion
            deleteFlag = True
            deleteSize = 1
        elif target_i is None and deleteFlag is True:
            #extending deletion
            deleteSize += 1
        elif target_i is not None and deleteFlag is True:
            #exiting deletion
            deleteFlag = False
            start = target_i - 1
            stop = target_i
            if mult3 is True and delete_size % 3 == 0:
                records.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,0,0", "A"))
            elif mult3 is False and delete_size % 3 == 0:
                records.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,0,0", "A"))
    if len(records) > 0:
        final.append(records)
    count += 1


    ##########
    #           0 1 2 3 4 5 6 7 8 9 1011121314151617
    # query     G T A T T - - T G G A C C T
    # target    G T A T T C T T G G A C C T A A G

aln = psl_lib.PslRow(simplePsl("+", 12, 0, 12, 290094216, 0, 12, [5,7], [0,5], [0,7], qName="query", tName="chr1"))
t = seq_lib.Transcript(bedLine("chr1", 0, 14, "query", 1000, "+", 0, 14, "128,0,0", 2, "5,7", "0,7").split())
a = seq_lib.Transcript(bedLine("chr1", 0, 12, "query", 1000, "+", 0, 12, "128,0,0", 1, 12, 0).split())



mult3 = False
records = []
final = []
count = 0
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    t = transcriptDict[aId]        
    records = []
    exonStarts = [x.start for x in a.exons]
    prevTargetPos = None
    for query_i in xrange(len(a)):
        if query_i in exonStarts:
            prevTargetPos = None
        target_i = aln.queryCoordinateToTarget(query_i)
        if target_i is None:
            #found deletion
            continue
        if prevTargetPos is not None and abs(target_i - prevTargetPos) != 1:
            #found insertion
            start = min(prevTargetPos, target_i) + 1
            stop = max(prevTargetPos, target_i)
            records.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,0,0", "A"))
        prevTargetPos = target_i
    if len(records) > 0:
        final.append(records)
    count += 1
