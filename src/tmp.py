import os

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from src.abstractClassifier import AbstractClassifier
a = AbstractClassifier(1,"datafiles/Rattus.filtered.psl",1,1,1,1,1,1,1,1)


transcripts = seq_lib.getTranscripts("datafiles/Rattus.gene-check.bed")
transcriptDict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
annotations = seq_lib.getTranscripts("datafiles/wgEncodeGencodeBasicVM2.gene-check.bed")
annotationDict = seq_lib.transcriptListToDict(annotations, noDuplicates=True)
alignments = psl_lib.readPsl("datafiles/Rattus.filtered.psl")
alignmentDict = psl_lib.getPslDict(alignments, noDuplicates=True)
seqDict = seq_lib.readTwoBit("datafiles/Rattus.2bit")


aId = "ENSMUST00000089106.2-1"

t = transcriptDict[aId]
e = t.exons[0]

t = transcriptDict[aId]
s = t.getCds(seqDict)[-3:]


valueDict = {}
for aId in alignmentDict.keys():
    if aId not in transcriptDict:
            continue
    t = transcriptDict[aId]
    if t.thickStart == t.thickStop == 0 or t.thickStop - t.thickStart < 3:
            continue
    s = t.getCds(seqDict)
    if not s.startswith("ATG"):
        valueDict[aId] = a.transcriptCoordinateToBed(t, 0, 2, "128,0,0", "BeginStart")
