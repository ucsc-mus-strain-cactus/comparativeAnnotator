import lib.sequence_lib as seq_lib


def insertionIterator(a, t, aln, mult3=False, inversion=False):
    """
    Target insertion:
    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA  

    Analyze a given annotation transcript, transcript and alignment for insertions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported. Set to None to report all.

    If inversion is True, inversions are ignored.
    """
    prev_target_i = None
    exon_starts = [x.start for x in a.exons]
    for query_i in xrange(len(a)):
        if query_i in exon_starts:
            # don't call a intron an insertion
            prev_target_i = None
        target_i = aln.queryCoordinateToTarget(query_i)
        if target_i is None:
            # deletion; ignore
            prev_target_i = target_i
            continue
        if prev_target_i is not None and abs(target_i - prev_target_i) != 1:
            # jumped over a insertion
            # make sure this is not an inversion
            if inversion is True:
                if t.chromosomeInterval.strand is True and target_i <= prev_target_i:
                    prev_target_i = target_i
                    continue
                if t.chromosomeInterval.strand is False and target_i >= prev_target_i:
                    prev_target_i = target_i
                    continue
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


def deletionIterator(a, t, aln, mult3=False, inversion=False):
    """
    Target deletion:
    query:   AATTATAAGCATGGA
    target:  AATTAT--GCATGGA

    Analyze a given annotation transcript, transcript and alignment for deletions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported.

    If inversion is True, inversions are ignored.
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
            # make sure this is not an inversion
            if inversion is True:
                if t.chromosomeInterval.strand is True and query_i <= prev_query_i:
                    prev_query_i = query_i
                    continue
                if t.chromosomeInterval.strand is False and query_i >= prev_query_i:
                    prev_query_i = query_i
                    continue
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

    """
    deletions = list(deletionIterator(a, t, aln, mult3=False, inversion=False))
    insertions = list(insertionIterator(a, t, aln, mult3=False, inversion=False))
    if len(deletions) == 0 and len(insertions) == 0:
        yield None
    for start, stop, span in sorted(deletions + insertions, key = lambda x: x[0]):
        yield start, stop, span


def codonPairIterator(a, t, aln, targetSeqDict, querySeqDict):
    """

    Inputs:
    Transcript objects representing the annotation (query) transcript and the target transcript.
    PslRow object that represents the alignment between the transcript objects.
    SeqDicts/TwoBitFileObjs that contain the genomic sequence for these two transcripts

    Yields matching codon pairs, taking into account indels. Out of frame pairs will not be returned

    Order is (target_cds_pos, target, query)
    
    """
    target_cds = t.getCds(targetSeqDict)
    query_cds = a.getCds(querySeqDict)
    if len(target_cds) == 0 or len(query_cds) == 0:
        yield None
    frame_shifts = list(frameShiftIterator(a, t, aln))
    if len(frame_shifts) > 0:
        frame_shifts = {start: size for start, stop, size in frame_shifts}
    frame_shift = False
    # iterate over the cds looking for codon pairs
    for target_cds_i in xrange(1, len(target_cds)):
        target_i = t.cdsCoordinateToChromosome(target_cds_i)
        query_i = aln.targetCoordinateToQuery(target_i)
        # the if statements below determine if we are moving in or out of frame
        if frame_shift is False and target_i in frame_shifts:
            frame_shift = True
            shift_size = frame_shifts[target_i]
        elif frame_shift is True and target_i in frame_shifts:
            shift_size += frame_shifts[target_i]
            if shift_size % 3 == 0:
                frame_shift = False
        # if we are in frame, we start yielding codon pairs
        if frame_shift is False and target_cds_i % 3 == 0:
            query_cds_i = a.transcriptCoordinateToCds(query_i)
            yield target_cds_i - 3, target_cds[target_cds_i - 3:target_cds_i], query_cds[query_cds_i - 3:query_cds_i]