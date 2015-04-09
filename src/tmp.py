import os

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from src.abstractClassifier import AbstractClassifier
from collections import defaultdict, Counter
from src.helperFunctions import *
from itertools import izip
from src.helperFunctions import *


transcripts = seq_lib.getTranscripts("pipeline_data/comparative/1411/transMap/results/geneCheck/C57B6NJ.gene-check.bed")
transcriptDict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
annotations = seq_lib.getTranscripts("pipeline_data/comparative/1411/transMap/data/wgEncodeGencodeBasicVM4.gene-check.bed")
annotationDict = seq_lib.transcriptListToDict(annotations, noDuplicates=True)
alignments = psl_lib.readPsl("pipeline_data/comparative/1411/transMap/results/filtered/C57B6NJ.filtered.psl")
alignmentDict = psl_lib.getPslDict(alignments, noDuplicates=True)
seqDict = seq_lib.readTwoBit("pipeline_data/comparative/1411/transMap/data/genomes/C57B6NJ.2bit")
refTwoBit = seq_lib.readTwoBit("pipeline_data/comparative/1411/transMap/data/genomes/C57B6J.2bit")
#refDict = seq_lib.getSequenceDict("../mouse_release_data/1411/C57B6J.fa")

aId = "ENSMUST00000121953.1-1" #disc1
a = annotationDict[aId[:-2]]
t = transcriptDict[aId]
aln = alignmentDict[aId]

aId = "ENSMUST00000112514.1-2" #start is out of frame
a = annotationDict[aId[:-2]]
t = transcriptDict[aId]
aln = alignmentDict[aId]


valueDict = {}
for aId, aln in alignmentDict.iteritems():
    if aln.strand == "+" and aln.tStart == 0 and aln.qStart != 0:
        valueDict[aId] = 1
    elif aln.strand == "-" and aln.tEnd == aln.tSize and aln.qEnd != aln.qSize:
        valueDict[aId] = 1
    else:
        valueDict[aId] = 0

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

aln = psl_lib.PslRow(psl_lib.simplePsl("+", 12, 0, 12, 290094216, 0, 12, [5,7], [0,5], [0,7], qName="query", tName="chr1"))
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


 ##########
    # DELETION
    #            0          11
    # ref        ATGATCCAATGA  query
    # exons       ****  ****
    # non ref    ATGATTAA--GA  target
    #            0          9
    #####

from lib.lib_tests import simplePsl
aln = simplePsl('+', 8, 0, 8, 10, 1, 9, [4, 1, 1], [0, 4, 7], [1, 7, 8], qName='ensmust0', tName='test_0_nr')
a = seq_lib.Transcript(['test_0_r', 1, 11, 'ensmust0', 0, '+', 1, 11, '128,0,0', 2, '4,4', '0,6'])
t = seq_lib.Transcript(['test_0_nr', 1, 9, 'ensmust0', 0, '+', 1, 9, '128,0,0', 2, '4,2', '0,6'])

    ##########
    # INSERTION
    #            0          9
    # ref        ATGATTAA--GA  query
    # exons       ****  ****
    # non ref    ATGATCCAATGA  target
    #            0          11
    #####

aln = simplePsl('+', 6, 0, 6, 12, 1, 11, [4, 1, 1], [0, 4, 5], [1, 7, 10], qName='ensmust0', tName='test_0_nr')
a = seq_lib.Transcript(['test_0_r', 1, 9, 'ensmust0', 0, '+', 1, 9, '128,0,0', 2, '4,2', '0,6'])
t = seq_lib.Transcript(['test_0_nr', 1, 11, 'ensmust0', 0, '+', 1, 11, '128,0,0', 2, '4,4', '0,6'])

    ##########
    # DELETION ON BOUNDARY
    #            0          11
    # ref        ATGATCCAATGA  query
    # exons       ****  ****
    # non ref    ATGATTA--AGA  target
    #            0          9
    #####

aln = simplePsl('+', 8, 0, 8, 10, 1, 9, [4, 2], [0, 6], [1, 7], qName='ensmust0', tName='test_0_nr')
a = seq_lib.Transcript(['test_0_r', 1, 11, 'ensmust0', 0, '+', 1, 11, '128,0,0', 2, '4,4', '0,6'])
t = seq_lib.Transcript(['test_0_nr', 1, 9, 'ensmust0', 0, '+', 1, 9, '128,0,0', 2, '4,2', '0,6'])



def getBed(t, rgb=None, name=None, start_offset=None, stop_offset=None):
    """
    Returns this transcript as a BED record with optional changes to rgb and name.
    If start_offset or stop_offset are set (chromosome coordinates), then this record will be changed to only 
    show results within that region, which is defined in chromosome coordinates.
    """
assert start_offset < stop_offset
if rgb is None:
    rgb = t.rgb
if name is not None:
    name += "/" + t.name
else:
    name = t.name
if start_offset is None and stop_offset is None:
    return [t.chrom, t.start, t.stop, name, t.score, convertStrand(t.strand), t.thickStart, t.thickStop, rgb, 
            t.blockCount, t.blockSizes, t.blockStarts]

def _moveStart(exonIntervals, blockCount, blockStarts, blockSizes, start, start_offset):
    toRemove = len([x for x in t.exonIntervals if x.start <= start_offset and x.stop <= start_offset])
    if toRemove > 0:
        blockCount -= toRemove
        blockSizes = blockSizes[toRemove:]
        start += blockStarts[toRemove]
        new_block_starts = [0]
        for i in xrange(toRemove, len(blockStarts) - 1):
            new_block_starts.append(blockStarts[i + 1] - blockStarts[i] + new_block_starts[-1])
        blockStarts = new_block_starts
    if start_offset > start:
        blockSizes[0] += start - start_offset
        blockStarts[1:] = [x + start - start_offset for x in blockStarts[1:]]
        start = start_offset
    return start, blockCount, blockStarts, blockSizes


def _moveStop(exonIntervals, blockCount, blockStarts, blockSizes, stop, start, stop_offset):
    toRemove = len([x for x in t.exonIntervals if x.stop >= stop_offset and x.start >= stop_offset])
    if toRemove > 0:
        blockCount -= toRemove
        blockSizes = blockSizes[:-toRemove]
        blockStarts = blockStarts[:-toRemove]
        stop = start + blockSizes[-1] + blockStarts[-1]
    if stop_offset < stop and stop_offset > start + blockStarts[-1]:
        blockSizes[-1] = stop_offset - start - blockStarts[-1] 
        stop = stop_offset
    return stop, blockCount, blockStarts, blockSizes
    

import os

import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
from src.abstractClassifier import AbstractClassifier
from collections import defaultdict, Counter
from itertools import izip


test = "chr1 100 10000 test 0 + 100 10000 128,0,0 5 100,150,100,1000,3000 0,1000,3000,4000,6900"
t = seq_lib.Transcript(test.split())
" ".join(map(str,t.getBed(start_offset=1150,name="start_1150")))
" ".join(map(str,t.getBed(stop_offset=7500,name="stop_7500")))
" ".join(map(str,t.getBed(stop_offset=4500,name="stop_4500")))
" ".join(map(str,t.getBed(stop_offset=4500,start_offset=1150,name="stop_4500_start_1150")))
" ".join(map(str,t.getBed(stop_offset=6500,start_offset=3000,name="stop_6500_start_3000")))
" ".join(map(str,t.getBed(stop_offset=6500,start_offset=3000,name="stop_4500_start_3000")))
" ".join(map(str,t.getBed(start_offset=1150, stop_offset=1200,name="start_1150_stop1200")))


blockCount = int(t.blockCount)
blockStarts = map(int, t.blockStarts.split(","))
blockSizes = map(int, t.blockSizes.split(","))
start = t.start
stop = t.stop
thickStart = t.thickStart
thickStop = t.thickStop

if start_offset is not None and start_offset > start:
    start, blockCount, blockStarts, blockSizes = _moveStart(t.exonIntervals, blockCount, blockStarts, blockSizes, start, start_offset)
if stop_offset is not None and stop_offset > stop and stop_offset > stop:
    stop, blockCount, blockStarts, blockSizes = _moveStop(t.exonIntervals, blockCount, blockStarts, blockSizes, stop, stop_offset)
if start > thickStart:
    thickStart = start
if stop < thickStop:
    thickStop = stop
return [t.chrom, start, stop, name, t.score, convertStrand(t.strand), thickStart, thickStop, rgb, blockCount,
        blockSizes, blockStarts]


valueDict = {}
for aId, aln in alignmentDict.iteritems():
    if aId not in transcriptDict:
        continue
    t = transcriptDict[aId]
    a = annotationDict[psl_lib.removeAlignmentNumber(aId)]
    s = list(frameShiftIterator(a, t, aln))
    if len(s) == 0:
        continue
    elif len(s) == 1:
        start, stop, size = s[0]
        valueDict[aId] = seq_lib.chromosomeCoordinateToBed(t, start, t.stop, "128,0,0", "A")
    else:
        tmp = []
        for i in xrange(1, len(s), 2):
            start = s[i-1][0]
            stop = s[i][1]
            tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, stop, "128,0,0", "A"))
        if i % 2 != 0:
            start = s[-1][0]
            tmp.append(seq_lib.chromosomeCoordinateToBed(t, start, t.stop, "128,0,0", "A"))
        valueDict[aId] = tmp




import os
import argparse

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger
from lib.general_lib import FileType, DirType, FullPaths, classesInModule
import lib.sqlite_lib as sql_lib

import src.classifiers, src.details, src.attributes
from src.constructDatabases import ConstructDatabases
from src.buildTracks import BuildTracks

# hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".fa"
gene_check_ext = ".gene-check.bed"

import src.classifiers, src.details, src.attributes

classifiers = classesInModule(src.classifiers)
details = classesInModule(src.details)
attributes = classesInModule(src.attributes)

def parseDir(genomes, targetDir, ext):
    pathDict = {}
    for g in genomes:
        path = os.path.join(targetDir, g + ext)
        if not os.path.exists(path):
            raise RuntimeError("{} does not exist".format(path))
        pathDict[g] = path
    return pathDict

geneCheckBedDict = parseDir(genomes_1411, "../mouse_release_data/1411", gene_check_ext)

p = BuildTracks("1411_output", genomes_1411, "AlignmentId", "../mouse_release_data/1411", geneCheckBedDict, "../mouse_release_data/wgEncodeGencodeBasicVM2.gene-check.bed")

p.run()


geneCheckBedDict = parseDir(genomes_1412, "../mouse_release_data/1412", gene_check_ext)

p2 = BuildTracks("1412_output", genomes_1412, "AlignmentId", "../mouse_release_data/1412", geneCheckBedDict, "../mouse_release_data/wgEncodeGencodeBasicVM2.gene-check.bed")

geneCheckBedDict = parseDir(genomes_test, "../mouse_release_data/1411", gene_check_ext)

p2.run()

p3 = BuildTracks("test_output", genomes_test, "AlignmentId", "../mouse_release_data/1411", geneCheckBedDict, "../mouse_release_data/wgEncodeGencodeBasicVM2.gene-check.bed")

bigBedDirs=`/bin/ls -1d 1412_output/bedfiles/* | paste -s -d ","`
python hal/assemblyHub/hal2assemblyHub.py /hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1412/cactus/1412.hal 1412_trackHub  --jobTree 1412_haljobtree --finalBigBedDirs ${bigBedDirs} --batchSystem=singleMachine --stats --shortLabel 1412 --longLabel 1412 --hub 1412 --maxThreads 20 &> 1412.log &

bigBedDirs=`/bin/ls -1d 1411_output/bedfiles/* | paste -s -d ","`
python hal/assemblyHub/hal2assemblyHub.py /cluster/home/jcarmstr/public_html/mouseBrowser_1411/1411.hal 1411_GPIP_trackHub --jobTree 1411_haljobtree --finalBigBedDirs ${bigBedDirs} --batchSystem=singleMachine --stats --shortLabel 1411_GPIP --longLabel 1411_GPIP --hub 1411_GPIP --maxThreads 20

export PYTHONPATH=./:${PYTHONPATH}
export PATH=./sonLib/bin:./submodules/jobTree/bin:./hal/bin/:${PATH}

bigBedDirs="test_output/bedfiles/transMap,test_output/bedfiles/everything,test_output/bedfiles/GPIP"
python hal/assemblyHub/hal2assemblyHub.py /cluster/home/jcarmstr/public_html/mouseBrowser_1411/1411.hal test_trackHub  --jobTree test_haljobtree --finalBigBedDirs ${bigBedDirs} --batchSystem=singleMachine --stats --shortLabel test --longLabel test --hub test --maxThreads 30 &> test.log &

for f in /cluster/home/ifiddes/ifiddes_hive/mus_strain_cactus/pipeline/results_1411/*/*details*; do n=`echo $f | cut -d "/" -f 9 | cut -d "." -f 2`; mkdir $n; bedToBigBed $f ~/ifiddes_hive/mouse_release_data/1411/$n.chrom.sizes $n/$n.bb; done


# testing rescuing starting frameshifts
# the below manual creation starts +1 bp
aln = psl_lib.PslRow("12 0 0 0 0 0 1 2 + chr1 290094216 4 20 chr1 290094216 4 20 3 2,1,9, 4,8,11, 4,8,11,")
a = seq_lib.Transcript("chr1 0 20 query 0 + 3 17 0 2 6,12, 0,8,".split())
t = seq_lib.Transcript("chr1 4 20 target 0 + 4 17 0 3 2,1,10, 0,4,6,".split())

#the below manual creation starts +2 bp
aln = psl_lib.PslRow("12 0 0 0 0 0 1 2 + chr1 290094216 5 20 chr1 290094216 5 20 3 1,1,10, 5,8,10, 5,8,10,")
a = seq_lib.Transcript("chr1 0 20 query 0 + 3 17 0 2 6,12, 0,8,".split())
t = seq_lib.Transcript("chr1 5 20 target 0 + 5 17 0 3 1,2,10, 0,2,5,".split())

#the below manual creation starts at +2 and has an insertion instead of deletion
aln = psl_lib.PslRow("13 0 0 0 1 2 0 0 + chr1 290094216 5 20 chr1 17 0 17 3 1,3,9, 5,8,11, 0,3,8,")
a = seq_lib.Transcript("chr1 0 20 query 0 + 3 17 0 2 6,12, 0,8,".split())
t = seq_lib.Transcript("chr1 0 17 target 0 + 0 14 0 1 17, 0,".split())
