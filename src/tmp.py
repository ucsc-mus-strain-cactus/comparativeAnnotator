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

#transcripts = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/C57B6NJ/transMapGencodeBasicVM4.gp")
transcripts = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/PAHARIEiJ/syn/transMapGencodeBasicVM4.gp")
transcriptDict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
annotations = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp")
annotationDict = seq_lib.transcriptListToDict(annotations, noDuplicates=True)
#alignments = psl_lib.readPsl("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/C57B6NJ/transMapGencodeBasicVM4.psl")
alignments = psl_lib.readPsl("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/transMap/2015-05-28/transMap/PAHARIEiJ/syn/transMapGencodeBasicVM4.psl")
alignmentDict = psl_lib.getPslDict(alignments, noDuplicates=True)
seqDict = seq_lib.getSequenceDict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1509/C57B6NJ.fa")
refSeqDict = seq_lib.getSequenceDict("/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1509/C57B6J.fa")
augustusTranscripts = seq_lib.getGenePredTranscripts("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1509/augustus/tmr/C57B6NJ.gp")
augustusTranscriptDict = seq_lib.transcriptListToDict(augustusTranscripts, noDuplicates=True)
