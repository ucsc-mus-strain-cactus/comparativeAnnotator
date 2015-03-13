import lib.sequence_lib as seq_lib


def insertionIterator(a, t, aln, mult3=False, inversion=False):
    """

    Target insertion:
    query:   AATTAT--GCATGGA
    target:  AATTATAAGCATGGA    

    Analyze a given annotation transcript, transcript and alignment for coding insertions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported. Set to None to report all.

    Inversion controls whether to filter out inversions or keep them (if you are concerned about frame shifts)

    """
    prev_target_i = None
    for query_i in xrange(len(a)):
        target_i = t.chromosomeCoordinateToTranscript(aln.queryCoordinateToTarget(query_i))
        if target_i is None:
            # deletion; ignore
            prev_target_i = target_i
            continue
        if prev_target_i is not None and abs(target_i - prev_target_i) != 1:
            # jumped over a deletion
            # make sure this is not an inversion
            if inversion is True:
                if t.chromosomeInterval.strand is True and target_i <= prev_target_i:
                    continue
                if t.chromosomeInterval.strand is False and target_i >= prev_target_i:
                    continue
            insertSize = abs(target_i - prev_target_i) - 1
            start = t.transcriptCoordinateToChromosome(target_i - 1)
            stop = t.transcriptCoordinateToChromosome(target_i)
            if t.chromosomeCoordinateToCds(start) is not None:
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

    Analyze a given annotation transcript, transcript and alignment for coding deletions.

    mult3 controls whether only multiple of 3 or only not multiple of 3 are reported.

    Inversion controls whether to filter out inversions or keep them (if you are concerned about frame shifts)

    """
    prev_query_i = None
    for target_i in xrange(len(t)):
        target_chrom_i = t.transcriptCoordinateToChromosome(target_i)
        query_i = aln.targetCoordinateToQuery(target_chrom_i)
        if query_i is None:
            # deletion; ignore
            prev_query_i = query_i
            continue
        if prev_query_i is not None and abs(query_i - prev_query_i) != 1:
            # jumped over a deletion
            # make sure this is not an inversion
            if inversion is True:
                if t.chromosomeInterval.strand is True and query_i <= prev_query_i:
                    continue
                if t.chromosomeInterval.strand is False and query_i >= prev_query_i:
                    continue            
            deleteSize = abs(query_i - prev_query_i) - 1
            start = stop = target_chrom_i
            if t.chromosomeCoordinateToCds(start) is not None:
                if mult3 is True and deleteSize % 3 == 0:
                    yield start, stop, deleteSize
                elif mult3 is False and deleteSize % 3 != 0:
                    yield start, stop, deleteSize
                elif mult3 is None:
                    yield start, stop, deleteSize
        prev_query_i = query_i


def firstIndel(a, t, aln, mult3=False, inversion=False):
    """

    Makes use of deletionIterator and insertionIterator to find the first indel in an alignment.

    If mult3 and inversion are False, this will be the first frameshift.

    """
    firstDeletion = next(deletionIterator(a, t, aln, mult3, inversion), None)
    firstInsertion = next(insertionIterator(a, t, aln, mult3, inversion), None)
    if firstDeletion is not None and firstInsertion is not None:
        # who happens first?
        if t.chromosomeInterval.strand is True:
            if firstInsertion > firstDeletion:
                return firstInsertion
            else:
                return firstDeletion                       
        else:
            if firstInsertion < firstDeletion:
                return firstInsertion
            else:
                return firstDeletion
    elif firstDeletion is not None:
        return firstDeletion
    elif firstInsertion is not None:                                       
        return firstInsertion
    return None


def codonPairIterator(a, t, aln, targetSeqDict, querySeqDict):
    """

    Inputs:
    Transcript objects representing the annotation (query) transcript and the target transcript.
    PslRow object that represents the alignment between the transcript objects.
    SeqDicts/TwoBitFileObjs that contain the genomic sequence for these two transcripts

    Yields matching codon pairs, taking into account indels. Out of frame pairs will not be returned

    Order is (target_pos, target, query)
    
    """
    target_cds = t.getCds(targetSeqDict)
    query_cds = a.getCds(querySeqDict)
    if len(target_cds) == 0 or len(query_cds) == 0:
        yield None
    # dicts mapping starts of indels to their size
    # deletions are stored as a negative value
    deletions = list(deletionIterator(a, t, aln, mult3=False, inversion=False))
    insertions = list(insertionIterator(a, t, aln, mult3=False, inversion=False))
    if len(deletions) > 0:
        deletions = {start : -int(size) for start, stop, size in deletions}
    if len(insertions) > 0:
        insertions = {start : int(size) for start, stop, size in insertions}
    frame_shift = False
    # iterate over the cds looking for codon pairs
    for target_cds_i in xrange(len(target_cds) - len(target_cds) % 3):
        target_i = t.cdsCoordinateToChromosome(target_cds_i)
        query_i = aln.targetCoordinateToQuery(target_i)
        # the if statements below determine if we are moving in or out of frame
        if frame_shift is False and target_i in insertions:
            # entering frame shift with insertion
            frame_shift = True
            shift_size = insertions[target_i]
        elif frame_shift is False and target_i in deletions:
            # entering frame shift with deletion
            frame_shift = True
            shift_size = deletions[target_i]
        elif frame_shift is True and target_i in insertions:
            # next frame shift - did we rescue?
            shift_size += insertions[target_i]
            if shift_size % 3 == 0:
                # rescued frame shift, we are back in frame
                frame_shift = False
                shift_size = 0
        elif frame_shift is True and target_i in deletions:
            # next frame shift - did we rescue?
            shift_size += deletions[target_i]
            if shift_size % 3 == 0:
                # rescued frame shift, we are back in frame
                frame_shift = False
                shift_size = 0
        # if we are in frame, we start yielding codon pairs
        if frame_shift is False and target_cds_i % 3 == 0:
            query_cds_i = a.transcriptCoordinateToCds(query_i)
            yield target_i, target_cds[target_cds_i:target_cds_i + 3], query_cds[query_cds_i:query_cds_i + 3]