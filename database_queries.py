"""
Uses the database schema to produce queries.
"""
from peewee import OperationalError
import time
from collections import Counter, defaultdict
from comparativeAnnotator.database_schema import ref_tables, tgt_tables, aug_tables, fetch_database


########################################################################################################################
## Initialize databases, helper functions
########################################################################################################################


def initialize_full_session(ref_genome, genome, db_path):
    """
    Returns all tables and establishes a session.
    """
    db = fetch_database()
    db.init(db_path)
    ref = ref_tables(ref_genome)
    tgt = tgt_tables(genome)
    aug = aug_tables(genome)
    return aug, tgt, ref


def initialize_tm_session(ref_genome, genome, db_path):
    """
    Returns all tables and establishes a session.
    """
    db = fetch_database()
    db.init(db_path)
    ref = ref_tables(ref_genome)
    tgt = tgt_tables(genome)
    return tgt, ref


def initialize_session(genome, db_path, table_fn):
    """
    Returns genome-specific tables with an established session.
    """
    db = fetch_database()
    db.init(db_path)
    return table_fn(genome)


def tgt_ref_join_aln_id(tgt, ref):
    """
    Produces a base joined table joining ref_attrs, tgt_attrs, ref_classify, tg_classify, returning IDs.
    """
    r = tgt.attrs.select(tgt.attrs.AlignmentId).\
        join(tgt.classify, on=(tgt.attrs.AlignmentId == tgt.classify.AlignmentId)).\
        join(ref.attrs, on=(ref.attrs.TranscriptId == tgt.classify.TranscriptId)).\
        join(ref.classify, on=(tgt.classify.TranscriptId == ref.classify.TranscriptId))
    return r


def aug_tgt_ref_join(aug, tgt, ref):
    """
    Produces a base joined table joining aug_attrs, aug_classify, ref_attrs, tgt_attrs, ref_classify, tg_classify,
    returning a combined object with all
    """
    r = aug.attrs.select(aug.attrs, aug.classify, tgt.attrs, tgt.classify).\
        join(aug.classify, on=(aug.attrs.AugustusAlignmentId == aug.classify.AugustusAlignmentId)).\
        join(tgt.classify, on=(aug.attrs.AlignmentId == tgt.classify.AlignmentId)).\
        join(tgt.attrs, on=(tgt.classify.AlignmentId == tgt.attrs.AlignmentId)).\
        join(ref.attrs, on=(ref.attrs.TranscriptId == tgt.classify.TranscriptId)).\
        join(ref.classify, on=(tgt.classify.TranscriptId == ref.classify.TranscriptId)).naive()
    return r


def tgt_ref_join(tgt, ref):
    """
    Produces a base joined table joining ref_attrs, tgt_attrs, ref_classify, tg_classify, returning a combined
    object with all
    """
    r = tgt.attrs.select(tgt.attrs, tgt.classify).\
        join(tgt.classify, on=(tgt.attrs.AlignmentId == tgt.classify.AlignmentId)).\
        join(ref.attrs, on=(ref.attrs.TranscriptId == tgt.classify.TranscriptId)).\
        join(ref.classify, on=(tgt.classify.TranscriptId == ref.classify.TranscriptId)).naive()
    return r


def add_biotype(r, ref, biotype):
    """
    Adds a biotype requirement to a query statement
    """
    return r.where(ref.attrs.GeneType == biotype)


def intron_inequality(r, tgt, ref):
    """
    Adds the original intron inequality to a select statement
    """
    return r.where(((2.0 * tgt.attrs.NumberMissingOriginalIntrons) <= (tgt.attrs.NumberIntrons - 1))
                   | (ref.attrs.NumberIntrons == 0))


def coding_classify(r, tgt, ref, passing):
    """
    Query that defines passing/excellent for coding genes, adding on to the noncoding classifiers.
    """
    r = r.where(((ref.classify.CdsUnknownSplice != 0) | (tgt.classify.CdsUnknownSplice == 0)),
                (tgt.classify.FrameShift == 0),
                (tgt.classify.CodingInsertions == 0),
                (tgt.classify.CodingDeletions == 0),
                (ref.classify.CdsGap != 0) | (tgt.classify.CdsGap == 0))

    if passing is False:
        r = r.where(((ref.classify.UtrGap != 0) | (tgt.classify.UtrGap == 0)),
                    ((ref.classify.BeginStart != 0) | (tgt.classify.BeginStart == 0)),
                    ((ref.classify.EndStop != 0) | (tgt.classify.EndStop == 0)),
                    ((ref.classify.InFrameStop != 0) | (tgt.classify.InFrameStop == 0)),
                    ((ref.classify.StartOutOfFrame != 0) | (tgt.classify.StartOutOfFrame == 0)),
                    ((ref.classify.ShortCds != 0) | (tgt.classify.ShortCds == 0)),
                    ((ref.classify.BadFrame != 0) | (tgt.classify.BadFrame == 0)),
                    (tgt.classify.HasOriginalStop == 0),
                    (tgt.classify.HasOriginalStart == 0),
                    (tgt.classify.UnknownGap == 0))
    return r


def noncoding_classify(r, tgt, ref, coverage, percent_unknown):
    """
    Constructs a query for noncoding classifiers. Adjust coverage and percent_unknown to adjust the amount of coverage
    and unknown bases allowed.
    """
    r = r.where(((ref.classify.UtrUnknownSplice != 0) | (tgt.classify.UtrUnknownSplice == 0)),
                (tgt.attrs.PercentUnknownBases <= percent_unknown),
                (tgt.attrs.AlignmentCoverage >= coverage))
    return r


def augustus_classify(r, aug, tgt, ref):
    """
    Constructs a query for augustus passing. Generally, allows Augustus to change things if things are wrong in
    the transMap or the reference to begin with. Also has a catch-all that allows any transcript with >95% coverage
    and >95% identity.
    TODO: don't hardcode the identity/coverage cutoffs.
    """
    # repeated requirements for both types of boundary movements
    boundaries = (tgt.attrs.AlignmentCoverage < 95.0) | (tgt.classify.CdsUnknownSplice > 0) | (tgt.classify.UtrUnknownSplice > 0)
    r = r.where((((aug.classify.NotSameStart == 0) | ((tgt.classify.HasOriginalStart != 0) |
                                                      (tgt.classify.StartOutOfFrame != 0) | (tgt.classify.BadFrame != 0) |
                                                      (ref.classify.BeginStart != 0))) &
                ((aug.classify.NotSameStop == 0) | ((tgt.classify.HasOriginalStop != 0) | (tgt.classify.BadFrame != 0) |
                                                    (ref.classify.EndStop != 0))) &
                ((aug.classify.NotSimilarTerminalExonBoundaries == 0) | (boundaries | (tgt.classify.UtrGap > 1))) &
                ((aug.classify.NotSimilarInternalExonBoundaries == 0) | (boundaries | (tgt.classify.CdsGap + tgt.classify.UtrGap > 2))) &
                ((aug.classify.ExonLoss < 2) & (aug.classify.AugustusParalogy == 0))) |
                ((aug.attrs.AugustusAlignmentCoverage >= 50.0) & (aug.attrs.AugustusAlignmentIdentity >= 95.0)))
    return r


def add_filter_chroms(r, ref, filter_chroms):
    """
    Take a list of chromosomes to filter out and add those to the query.
    """
    for c in filter_chroms:
        r = r.where(ref.attrs.SourceChrom != c)
    return r


def add_best_cov(r, tgt):
    """
    Adds the bestcoverage requirement to a query.
    """
    return r.where(tgt.attrs.HighestCovAln == True)


########################################################################################################################
## Queries
########################################################################################################################


def execute_query(query, timeout=6000, interval=10):
    """
    Handle exceptions to be more informative, as well as handles database timeouts. Will try again every interval
    """
    start_time = time.time()
    while time.time() - start_time <= timeout:
        try:
            return list(query.execute())
        except OperationalError, e:
            if 'locked' in e:
                time.sleep(interval)
            else:
                raise OperationalError('Error. Original message: {}'.format(e))
    raise OperationalError('Error: database still locked after {} seconds'.format(timeout))


def get_fail_pass_excel_ids(ref_genome, genome, db_path, biotype=None, filter_chroms=None, best_cov_only=False):
    all_ids = get_aln_ids(ref_genome, genome, db_path, biotype, best_cov_only)
    pass_ids = trans_map_eval(ref_genome, genome, db_path, biotype=biotype, filter_chroms=filter_chroms, passing=True,
                              best_cov_only=best_cov_only)
    excel_ids = trans_map_eval(ref_genome, genome, db_path, biotype=biotype, filter_chroms=filter_chroms, passing=False,
                               best_cov_only=best_cov_only)
    pass_specific_ids = pass_ids - excel_ids
    fail_ids = all_ids - pass_ids
    assert len(fail_ids) + len(pass_specific_ids) + len(excel_ids) == len(all_ids)
    return excel_ids, pass_specific_ids, fail_ids


def trans_map_eval(ref_genome, genome, db_path, biotype=None, filter_chroms=None, passing=False, best_cov_only=False):
    """
    Returns the set of AlignmentIds that pass the classifiers.
    """
    tgt, ref = initialize_tm_session(ref_genome, genome, db_path)
    r = tgt_ref_join_aln_id(tgt, ref)
    r = intron_inequality(r, tgt, ref)
    if passing is True:
        r = noncoding_classify(r, tgt, ref, coverage=90.0, percent_unknown=5.0)
    else:
        r = noncoding_classify(r, tgt, ref, coverage=99.0, percent_unknown=1.0)
    if biotype == 'protein_coding':
        r = coding_classify(r, tgt, ref, passing)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if filter_chroms is not None:
        r = add_filter_chroms(r, ref, filter_chroms)
    if best_cov_only is True:
        r = add_best_cov(r, tgt)
    return set([x[0] for x in execute_query(r.tuples())])


def augustus_eval(ref_genome, genome, db_path, biotype=None, filter_chroms=None):
    """
    Returns the set of AlignmentIds that pass the classifiers.
    """
    aug, tgt, ref = initialize_full_session(ref_genome, genome, db_path)
    r = aug_tgt_ref_join(aug, tgt, ref)
    r = augustus_classify(r, aug, tgt, ref)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if filter_chroms is not None:
        r = add_filter_chroms(r, ref, filter_chroms)
    return set([x[0] for x in execute_query(r.tuples())])


def get_aln_ids(ref_genome, genome, db_path, biotype=None, best_cov_only=False):
    """
    Gets the alignment IDs in the database. Filters on biotype if set.
    """
    tgt, ref = initialize_tm_session(ref_genome, genome, db_path)
    r = tgt_ref_join_aln_id(tgt, ref)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if best_cov_only is True:
        r = add_best_cov(r, tgt)
    return set([x[0] for x in execute_query(r.tuples())])


def get_ref_ids(ref_genome, db_path, biotype=None):
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    return set([x[0] for x in execute_query(r.tuples())])


def get_column(genome, ref_genome, db_path, col, biotype=None, best_cov_only=True):
    tgt, ref = initialize_tm_session(ref_genome, genome, db_path)
    r = tgt_ref_join(tgt, ref)
    r = r.select(eval(col))
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if best_cov_only is True:
        r = add_best_cov(r, tgt)
    return [x[0] for x in execute_query(r.tuples())]


def paralogy(genome, db_path, biotype=None):
    """
    Returns a counter of paralogy.
    """
    tgt = initialize_session(genome, db_path, tgt_tables)
    r = tgt.attrs.select(tgt.attrs.Paralogy)
    return [x[0] + 1 for x in execute_query(r.tuples())]


def get_transcript_gene_map(ref_genome, db_path, biotype=None):
    """
    Returns a dict mapping transcript names to gene names.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId, ref.attrs.GeneId)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    result = {}
    for tx_id, gene_id in execute_query(r.tuples()):
        result[tx_id] = gene_id
    return result


def get_transcript_biotype_map(ref_genome, db_path):
    """
    Returns a dictionary mapping transcript IDs to transcript biotypes.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId, ref.attrs.TranscriptType)
    result = {}
    for tx_id, biotype in execute_query(r.tuples()):
        result[tx_id] = biotype
    return result


def get_gene_transcript_map(ref_genome, db_path, biotype=None):
    """
    Returns a dict mapping transcript names to gene names.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId, ref.attrs.GeneId)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    result = defaultdict(list)
    for tx_id, gene_id in execute_query(r.tuples()):
        result[gene_id].append(tx_id)
    return result


def get_biotypes(ref_genome, db_path):
    """
    Returns the set of biotypes present in this reference set.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptType)
    return set(x[0] for x in execute_query(r.tuples()))


def get_rows(ref_genome, genome, db_path, mode='transMap', biotype=None):
    """
    Returns a combined iterable of all rows in the classify-ref table set.
    """
    if mode == 'transMap':
        tgt, ref = initialize_tm_session(ref_genome, genome, db_path)
        r = tgt_ref_join(tgt, ref)
    elif mode == 'AugustusTMR' or mode == 'AugustusTM':
        aug, tgt, ref = initialize_full_session(ref_genome, genome, db_path)
        r = aug_tgt_ref_join(aug, tgt, ref)
    else:
        raise Exception("bad programmer")
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    try:
        return execute_query(r.naive())
    except Exception, e:
        assert False, (mode, r, e)


def get_row_dict(ref_genome, genome, db_path, mode, biotype=None):
    """
    Wraps get_rows, returning a dictionary mapping the unique ID to the row.
    """
    col = 'x.AlignmentId' if mode == 'transMap' else 'x.AugustusAlignmentId'
    return {eval(col): x for x in get_rows(ref_genome, genome, db_path, mode, biotype)}
