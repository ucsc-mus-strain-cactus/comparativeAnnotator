"""
Uses the database schema to produce queries.
"""
from collections import Counter, defaultdict
from comparativeAnnotator.database_schema import ref_tables, tgt_tables, fetch_database


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


def tgt_ref_join(tgt, ref):
    """
    Produces a base joined table joining ref_attrs, tgt_attrs, ref_classify, tg_classify, returning a combined
    object with all
    """
    r = tgt.attrs.select(tgt.attrs, tgt.classify, ref.attrs, ref.classify).\
        join(tgt.classify, on=(tgt.attrs.AlignmentId == tgt.classify.AlignmentId)).\
        join(ref.attrs, on=(ref.attrs.TranscriptId == tgt.classify.TranscriptId)).\
        join(ref.classify, on=(tgt.classify.TranscriptId == ref.classify.TranscriptId)).naive()
    return r


def add_biotype(r, ref, biotype):
    """
    Adds a biotype requirement to a query statement
    """
    return r.where(ref.attrs.TranscriptType == biotype)


def intron_inequality(r, tgt, ref):
    """
    Adds the missing original intron inequality to a select statement
    """
    return r.where(((2.0 * tgt.attrs.NumberMissingOriginalIntrons) <= (tgt.attrs.NumberIntrons - 1))
                    | (ref.attrs.NumberIntrons == 0))


def coding_classify(r, tgt, ref, passing):
    r = r.where(((ref.classify.CdsUnknownSplice != 0) | (tgt.classify.CdsUnknownSplice == 0)),
                (tgt.classify.FrameShift == 0),
                (tgt.classify.CodingInsertions == 0),
                tgt.classify.CodingDeletions == 0)
    if passing is False:
        r = r.where(((ref.classify.CdsGap != 0) | (ref.classify.CdsGap == 0)),
                    ((ref.classify.BeginStart != 0) | (ref.classify.BeginStart == 0)),
                    ((ref.classify.EndStop != 0) | (ref.classify.EndStop == 0)),
                    ((ref.classify.InFrameStop != 0) | (ref.classify.InFrameStop == 0)),
                    ((ref.classify.StartOutOfFrame != 0) | (ref.classify.StartOutOfFrame == 0)),
                    ((ref.classify.ShortCds != 0) | (ref.classify.ShortCds == 0)),
                    ((ref.classify.BadFrame != 0) | (ref.classify.BadFrame == 0)),
                    (tgt.classify.HasOriginalStop == 0),
                    (tgt.classify.HasOriginalStart == 0))
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


def add_filter_chroms(r, ref, filter_chroms):
    """
    Take a list of chromosomes to filter out and add those to the query.
    """
    for c in filter_chroms:
        r = r.where(ref.attrs.RefChrom != c)
    return r


def add_best_cov(r, tgt):
    """
    Adds the bestcoverage requirement to a query.
    """
    return r.where(tgt.attrs.HighestCovAln == True)


########################################################################################################################
## Queries
########################################################################################################################


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
    tgt, ref = initialize_full_session(ref_genome, genome, db_path)
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
    return set([x[0] for x in r.tuples().execute()])


def get_aln_ids(ref_genome, genome, db_path, biotype=None, best_cov_only=False):
    """
    Gets the alignment IDs in the database. Filters on biotype if set.
    """
    tgt, ref = initialize_full_session(ref_genome, genome, db_path)
    r = tgt_ref_join_aln_id(tgt, ref)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if best_cov_only is True:
        r = add_best_cov(r, tgt)
    return set([x[0] for x in r.tuples().execute()])


def get_ref_ids(ref_genome, db_path, biotype=None):
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    return set([x[0] for x in r.tuples().execute()])


def get_column(genome, ref_genome, db_path, col, biotype=None, best_cov_only=True):
    tgt, ref = initialize_full_session(ref_genome, genome, db_path)
    r = tgt_ref_join(tgt, ref)
    r = r.select(eval(col))
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    if best_cov_only is True:
        r = add_best_cov(r, tgt)
    return [x[0] for x in r.tuples().execute()]


def paralogy(genome, db_path, biotype=None):
    """
    Returns a counter of paralogy.
    """
    tgt = initialize_session(genome, db_path, tgt_tables)
    r = tgt.attrs.select(tgt.attrs.Paralogy).tuples().execute()
    return [x[0] + 1 for x in r]


def get_transcript_gene_map(ref_genome, db_path, biotype=None):
    """
    Returns a dict mapping transcript names to gene names.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptId, ref.attrs.GeneId)
    if biotype is not None:
        r = add_biotype(r, ref, biotype)
    result = {}
    for tx_id, gene_id in r.tuples().execute():
        result[tx_id] = gene_id
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
    for tx_id, gene_id in r.tuples().execute():
        result[gene_id].append(tx_id)
    return result


def get_biotypes(ref_genome, db_path):
    """
    Returns the set of biotypes present in this reference set.
    """
    ref = initialize_session(ref_genome, db_path, ref_tables)
    r = ref.attrs.select(ref.attrs.TranscriptType)
    return set(x[0] for x in r.tuples().execute())

# this is how peewee can select items from multiple databases at once
#q=a.select(a, c).join(c, on=(a.TranscriptId == c.TranscriptId)).naive().execute()
