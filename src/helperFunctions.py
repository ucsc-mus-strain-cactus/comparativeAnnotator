from collections import defaultdict
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

def insertionIterator(a, t, aln, mult3=None):
    """
    Target insertion:
    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA  

    Analyze a given annotation transcript, transcript and alignment for insertions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported. Set to None to report all.
    """
    prev_target_i = None
    exon_starts = [x.start for x in a.exons]
    for query_i in xrange(len(a)):
        if query_i in exon_starts:
            # don't call a intron an insertion
            prev_target_i = None
            continue
        target_i = aln.queryCoordinateToTarget(query_i)
        if target_i is None:
            # deletion; ignore
            prev_target_i = target_i
            continue
        if prev_target_i is not None and abs(target_i - prev_target_i) != 1:
            # jumped over a insertion
            insertSize = abs(target_i - prev_target_i) - 1
            start = min(prev_target_i, target_i) + 1
            stop = max(prev_target_i, target_i)
            if mult3 is True and insertSize % 3 == 0:
                yield start, stop, insertSize
            elif mult3 is False and insertSize % 3 != 0:
                yield start, stop, insertSize
            elif mult3 is None:
                yield start, stop, insertSize
        prev_target_i = target_i


def deletionIterator(a, t, aln, mult3=None):
    """
    Target deletion:
    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA

    Analyze a given annotation transcript, transcript and alignment for deletions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported.
    """
    prev_query_i = None
    for target_i in xrange(len(t)):
        target_chrom_i = t.transcriptCoordinateToChromosome(target_i)
        query_i = aln.targetCoordinateToQuery(target_chrom_i)
        if query_i is None:
            # insertion; ignore
            prev_query_i = query_i
            continue
        if prev_query_i is not None and abs(query_i - prev_query_i) != 1:
            # jumped over a deletion
            deleteSize = abs(query_i - prev_query_i) - 1
            start = stop = target_chrom_i
            if mult3 is True and deleteSize % 3 == 0:
                yield start, stop, -deleteSize
            elif mult3 is False and deleteSize % 3 != 0:
                yield start, stop, -deleteSize
            elif mult3 is None:
                yield start, stop, -deleteSize
        prev_query_i = query_i


def frameShiftIterator(a, t, aln):
    """
    Yields frameshift-causing mutations. These are defined as non mult3 indels within CDS.
    These are returned in CDS coordinates.
    """
    initial = []
    start_shift = - (a.transcriptCoordinateToCds(aln.targetCoordinateToQuery(t.cdsCoordinateToChromosome(0))) % 3)
    if start_shift != 0:
        initial.append([0, start_shift])
    deletions = [[t.chromosomeCoordinateToCds(x), z] for x, y, z in deletionIterator(a, t, aln, mult3=False) if x > t.thickStart and y < t.thickStop]
    if t.strand is True:
        insertions = [[t.chromosomeCoordinateToCds(x - 1), z] for x, y, z in insertionIterator(a, t, aln, mult3=False) if x > t.thickStart and y < t.thickStop]
        d = defaultdict(int)
        for p, s in deletions + insertions + initial:
            d[p] += s
    else:
        insertions = [[t.chromosomeCoordinateToCds(y), z] for x, y, z in insertionIterator(a, t, aln, mult3=False) if x > t.thickStart and y < t.thickStop]
        d = defaultdict(int)
        for p, s in deletions + insertions + initial:
            d[p] += s
    combined = sorted(d.iteritems(), key = lambda x: x[0])
    for cds_pos, span in combined:
        yield cds_pos, span


def codonPairIterator(a, t, aln, targetSeqDict, querySeqDict):
    """
    Inputs:
    Transcript objects representing the annotation (query) transcript and the target transcript.
    PslRow object that represents the alignment between the transcript objects.
    SeqDicts/TwoBitFileObjs that contain the genomic sequence for these two transcripts

    Order is (target_cds_pos, target, query)
    """
    target_cds = t.getCds(targetSeqDict)
    query_cds = a.getCds(querySeqDict)
    for i in xrange(0, a.getCdsLength(), 3):
        target_cds_positions = [t.chromosomeCoordinateToCds(aln.queryCoordinateToTarget(a.cdsCoordinateToTranscript(j))) 
                                for j in xrange(i, i + 3)]
        if None in target_cds_positions:
            continue
        assert all([target_cds_positions[2] - target_cds_positions[1] == 1, target_cds_positions[1] - target_cds_positions[0] == 1, target_cds_positions[2] - target_cds_positions[0] == 2])
        target_codon = target_cds[target_cds_positions[0]:target_cds_positions[0] + 3]
        query_codon = query_cds[i:i + 3]
        yield target_cds_positions[0], target_codon, query_codon