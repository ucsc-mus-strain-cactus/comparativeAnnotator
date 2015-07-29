import os

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from src.abstractClassifier import AbstractClassifier
from collections import defaultdict, Counter
from src.helperFunctions import *
from src.classifiers import *
from itertools import izip
import intervaltree

transcripts = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/NODShiLtJ/transMapGencodeBasicVM4.gp")
transcriptDict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
annotations = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp")
annotationDict = seq_lib.transcriptListToDict(annotations, noDuplicates=True)
alignments = psl_lib.readPsl("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/NODShiLtJ/transMapGencodeBasicVM4.psl")
alignmentDict = psl_lib.getPslDict(alignments, noDuplicates=True)
seqDict = seq_lib.getSequenceDict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/NODShiLtJ.fa")
refSeqDict = seq_lib.getSequenceDict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/C57B6J.fa")
augustusTranscripts = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/augustus/tmr/NODShiLtJ.gp")
augustusTranscriptDict = seq_lib.transcriptListToDict(augustusTranscripts, noDuplicates=True)

aId = "ENSMUST00000068580.3-1"
aug_aId = "aug-ENSMUST00000068580.3-1"
aug_t = augustusTranscriptDict[aug_aId]
t = transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]


aId = "ENSMUST00000072079.8-4"
aug_aId = "aug-ENSMUST00000072079.8-4"
aug_t = augustusTranscriptDict[aug_aId]
t = transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]

t = Nonsynonymous("NODShiLtJ", "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/NODShiLtJ/transMapGencodeCompVM4.psl",
    "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/assemblies/1504/NODShiLtJ.fa",
    "/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/C57B6J.fa",
    "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeCompVM4.gp",
    "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv",
    "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/NODShiLtJ/transMapGencodeCompVM4.gp",
    "C57B6J", "AlignmentId", "./tmp")


aug_t_intervals = [seq_lib.ChromosomeInterval('1', 0, 10, True), seq_lib.ChromosomeInterval('1', 20, 30, True)]
merged_t_intervals = [seq_lib.ChromosomeInterval('1', 0, 8, True), seq_lib.ChromosomeInterval('1', 15, 18, True), seq_lib.ChromosomeInterval('1', 20, 30, True)]

shortIntronSize=30
wiggleRoom=10
classifyDict = {}
detailsDict = defaultdict(list)
strand = []
for aug_aId, aug_t in augustusTranscriptDict.iteritems():
    t = transcriptDict[psl_lib.removeAugustusAlignmentNumber(aug_aId)]
    if aug_t.strand != t.strand:
        strand.append(aug_aId)
        continue
    merged_t_intervals = seq_lib.gapMergeIntervals(t.exonIntervals, gap=shortIntronSize)
    aug_t_intervals = aug_t.exonIntervals
    for interval in merged_t_intervals:
        if seq_lib.intervalNotWithinWiggleRoomIntervals(aug_t_intervals, interval, wiggleRoom):
            classifyDict[aug_aId] = 1
            detailsDict[aug_aId].append(interval.getBed("A", "A"))
    if aug_aId not in classifyDict:
        classifyDict[aug_aId] = 0
