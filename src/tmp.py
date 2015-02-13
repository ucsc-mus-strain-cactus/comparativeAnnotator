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


insertFlag = False
startFlag = False
mult3 = False
records = []
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    for exon in annotationDict[psl_lib.removeAlignmentNumber(aId)].exons:
        for i in xrange(exon.start, exon.stop):
            #make sure we have actually entered the alignment before we call insertions
            if startFlag is False and aln.queryCoordinateToTarget(i) == None:
                continue
            # we have found an insertion
            if insertFlag is False and aln.queryCoordinateToTarget(i) == None:
                insertSize = 1
                insertFlag = True
            #insertion continues
            elif insertFlag is True and aln.queryCoordinateToTarget(i) == None:
                insertSize += 1
            #exiting insertion
            elif insertFlag is True and aln.queryCoordinateToTarget(i) != None:
                if (a.transcriptCoordinateToCds(i) is not None or
                            a.transcriptCoordinateToCds(i + insertSize) is not None):
                    t = transcriptDict[aId]
                    insertFlag = False
                    start = aln.queryCoordinateToTarget(i - insertSize - 1)
                    stop = aln.queryCoordinateToTarget(i)                    
                    if insertSize % 3 == 0 and mult3 == True:
                        if t.chromosomeInterval.strand is not True:
                            start, stop = stop, start
                        records.append(
                            seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,128,128", "A"))
                    elif insertSize % 3 != 0 and mult3 == False:
                        if t.chromosomeInterval.strand is not True:
                            start, stop = stop, start
                        records.append(
                            seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,128,128", "A"))
                insertSize = 0
            else:
                startFlag = True
                