"""
Convenience library for interfacing with a sqlite database.
"""
import os
import sys
import lib.psl_lib as psl_lib
import lib.general_lib as general_lib
from collections import defaultdict
import sqlite3 as sql
import pandas as pd
import etc.config

__author__ = "Ian Fiddes"


class ExclusiveSqlConnection(object):
    """meant to be used with a with statement to ensure proper closure"""
    def __init__(self, path, timeout=1200):
        self.path = path
        self.timeout = timeout

    def __enter__(self):
        self.con = sql.connect(self.path, timeout=self.timeout, isolation_level="EXCLUSIVE")
        try:
            self.con.execute("BEGIN EXCLUSIVE")
        except sql.OperationalError:
            raise RuntimeError("Database still locked after {} seconds.".format(self.timeout))
        return self.con

    def __exit__(self, exception_type, exception_val, trace):
        self.con.commit()
        self.con.close()


def attach_database(con, path, name):
    """
    Attaches another database found at path to the name given in the given connection.
    """
    con.execute("ATTACH DATABASE '{}' AS {}".format(path, name))


def open_database(path, timeout=1200):
    con = sql.connect(path, timeout=timeout)
    cur = con.cursor()
    return con, cur


def attach_databases(comp_ann_path, mode):
    """
    Attaches all of the databases needed for this execution mode.
    """
    assert mode in ["reference", "augustus", "transMap"]
    classify_path = os.path.join(comp_ann_path, "classify.db")
    con, cur = open_database(classify_path)
    details_path = os.path.join(comp_ann_path, "details.db")
    attach_database(con, details_path, "details")
    attr_path = os.path.join(comp_ann_path, "attributes.db")
    attach_database(con, attr_path, "attributes")
    if mode == "augustus":
        aug_classify_path = os.path.join(comp_ann_path, "augustus_classify.db")
        aug_details_path = os.path.join(comp_ann_path, "augustus_details.db")
        aug_attributes_path = os.path.join(comp_ann_path, "augustus_attributes.db")
        attach_database(con, aug_classify_path, "augustus")
        attach_database(con, aug_details_path, "augustus_details")
        attach_database(con, aug_attributes_path, "augustus_attributes")
    return con, cur


def load_data(con, genome, columns, primary_key="AlignmentId", table="main"):
    """
    Use pandas to load a sql query into a dataframe.
    """
    columns = ",".join(columns)
    query = "SELECT {},{} FROM {}.'{}'".format(primary_key, columns, table, genome)
    return pd.read_sql_query(query, con, index_col=primary_key)


def get_query_ids(cur, query):
    """
    Returns a set of aln_ids which are OK based on the definition of OK in config.py that made this query.
    In other words, expects a query of the form SELECT AlignmentId FROM stuff
    """
    try:
        return {x[0] for x in cur.execute(query)}
    except sql.OperationalError, exc:
        raise RuntimeError("query failed: {}.\nOriginal error message: {}".format(query, exc))


def get_query_dict(cur, query):
    """
    Returns a set of aln_ids which are OK based on the definition of OK in config.py that made this query.
    In other words, expects a query of the form SELECT AlignmentId,<other stuff> FROM stuff
    """
    try:
        return {x[0]: x[1] if len(x) == 2 else x[1:] for x in cur.execute(query)}
    except sql.OperationalError, exc:
        raise RuntimeError("query failed: {}.\nOriginal error message: {}".format(query, exc))


def get_non_unique_query_dict(cur, query, flatten=True):
    """
    Same as get_query_dict, but has no guarantee of uniqueness. Therefore, the returned data structure is a
    defaultdict(list)
    """
    d = defaultdict(list)
    try:
        for r in cur.execute(query):
            d[r[0]].append(r[1:])
    except sql.OperationalError, exc:
        raise RuntimeError("query failed: {}.\nOriginal error message: {}".format(query, exc))
    return general_lib.flatten_defaultdict_list(d) if flatten is True else d


def get_multi_index_query_dict(cur, query, num_indices=2):
    """
    Similar to get_non_unique_query_dict, but instead produces a layered dict of dicts where the number of layers
    is the num_indicies. Expects the query to be in the form Index1,Index2,IndexN,data and will produce the data
    structure {Index1: {Index2: data}}.
    """
    d = {}
    try:
        for r in cur.execute(query):
            last_layer = d
            for i in xrange(num_indices - 1):
                if r[i] not in last_layer:
                    last_layer[r[i]] = {}
                last_layer = last_layer[r[i]]
            last_layer[r[i + 1]] = r[i + 2:]
    except sql.OperationalError, exc:
        raise RuntimeError("query failed: {}.\nOriginal error message: {}".format(query, exc))
    return d


def get_biotype_aln_ids(cur, genome, biotype):
    """
    Returns a set of aln_ids which are members of the biotype.
    """
    query = "SELECT AlignmentId FROM attributes.'{0}' WHERE TranscriptType = '{1}' AND GeneType = '{1}'"
    query = query.format(genome, biotype)
    return get_query_ids(cur, query)


def get_biotype_ids(cur, ref_genome, biotype, mode="Transcript", filter_chroms=None):
    """
    Returns a set of transcript_ids which are members of the biotype. Category should be Transcript or Gene.
    In gene mode returns the gene_id instead of transcript_id.
    """
    assert mode in ["Transcript", "Gene"]
    query = ("SELECT {0}Id FROM main.'{1}' JOIN attributes.'{1}' USING (TranscriptId) WHERE "
             "TranscriptType = '{2}' AND GeneType = '{2}'")
    if filter_chroms is not None:
        for filter_chrom in filter_chroms:
            query += " AND refChrom != '{}'".format(filter_chrom)
    query = query.format(mode, ref_genome, biotype)
    return get_query_ids(cur, query)


def get_all_biotypes(cur, ref_genome, gene_level=False):
    """
    Returns all biotypes in the attribute database.
    """
    r = "Gene" if gene_level else "Transcript"
    query = "SELECT {}Type FROM attributes.'{}'".format(r, ref_genome)
    return get_query_ids(cur, query)


def get_ids_by_chromosome(cur, genome, chromosomes=("Y", "chrY")):
    """
    Returns all AlignmentIds in the attribute database whose source chromosome is in chromomsomes
    """
    base_query = "SELECT AlignmentId FROM attributes.'{}' WHERE {}"
    filt = ["SourceChrom = '{}'".format(x) for x in chromosomes]
    filt = " OR ".join(filt)
    query = base_query.format(genome, filt)
    return get_query_ids(cur, query)


def get_transcript_gene_map(cur, ref_genome, biotype=None, filter_chroms=None):
    """
    Returns a dictionary mapping transcript IDs to gene IDs
    """
    if biotype is None:
        query = "SELECT transcriptId,geneId FROM attributes.'{0}'"
    else:
        query = "SELECT transcriptId,geneId FROM attributes.'{0}' WHERE transcriptType = '{1}' AND geneType = '{1}'"
    if filter_chroms is not None:
        add = "AND".join([" refChrom != '{}' ".format(filter_chrom) for filter_chrom in filter_chroms])
        if biotype is None:
            query += " WHERE "
        else:
            query += " AND "
        query += add
    query = query.format(ref_genome, biotype)
    return get_query_dict(cur, query)


def get_gene_transcript_map(cur, ref_genome, biotype=None, filter_chroms=None):
    """
    Returns a dictionary mapping all gene IDs to their respective transcripts
    """
    if biotype is None:
        query = "SELECT geneId,transcriptId FROM attributes.'{0}'"
    else:
        query = "SELECT geneId,transcriptId FROM attributes.'{0}' WHERE transcriptType = '{1}' AND geneType = '{1}'"
    if filter_chroms is not None:
        add = "AND".join([" refChrom != '{}' ".format(filter_chrom) for filter_chrom in filter_chroms])
        if biotype is None:
            query += " WHERE "
        else:
            query += " AND "
        query += add
    query = query.format(ref_genome, biotype)
    return get_non_unique_query_dict(cur, query)


def get_gene_biotype_map(cur, ref_genome):
    """
    Returns a dictionary mapping all gene IDs to their respective biotypes
    """
    query = "SELECT geneId,geneType FROM attributes.'{}'".format(ref_genome)
    return get_query_dict(cur, query)


def get_transcript_biotype_map(cur, ref_genome):
    """
    Returns a dictionary mapping all transcript IDs to their respective biotypes
    """
    query = "SELECT transcriptId,transcriptType FROM attributes.'{}'".format(ref_genome)
    return get_query_dict(cur, query)


def get_fail_good_pass_ids(cur, ref_genome, genome, biotype, best_cov_only=True):
    """
    Returns the IDs categorized as fail, good_specific, pass
    """
    good_query = etc.config.transMapEval(ref_genome, genome, biotype, good=True)
    pass_query = etc.config.transMapEval(ref_genome, genome, biotype, good=False)
    if best_cov_only is True:
        best_covs = highest_cov_aln(cur, genome)
        best_ids = set(zip(*best_covs.itervalues())[0])
    else:
        class Universe:
            def __contains__(_,x): return True
        best_ids = Universe()
    good_ids = {x for x in get_query_ids(cur, good_query) if x in best_ids}
    pass_ids = {x for x in get_query_ids(cur, pass_query) if x in best_ids}
    all_ids = {x for x in get_biotype_aln_ids(cur, genome, biotype) if x in best_ids}
    fail_ids = all_ids - good_ids
    good_specific_ids = good_ids - pass_ids
    return fail_ids, good_specific_ids, pass_ids


def write_dict(data_dict, database_path, table, index_label="AlignmentId"):
    """
    Writes a dict of dicts to a sqlite database.
    """
    df = pd.DataFrame.from_dict(data_dict)
    df = df.sort_index()
    with ExclusiveSqlConnection(database_path) as con:
        df.to_sql(table, con, if_exists="replace", index_label=index_label)


def write_csv(csv_path, database_path, table, sep=",", index_col=0, header=0, index_label="AlignmentId"):
    """
    Writes a csv/tsv file to a sqlite database. Assumes that this table has a header
    """
    df = pd.read_table(csv_path, sep=sep, index_col=index_col, header=header)
    df = df.sort_index()
    with ExclusiveSqlConnection(database_path) as con:
        df.to_sql(table, con, if_exists="replace", index_label=index_label)


def collapse_details_dict(details_dict):
    """
    Collapses a details dict into a string so that the database record is a string. Properly accounts for the three
    possibilities: the record is a single BED entry (a string), the record is a list of BED entries, or the record
    is empty but has been accessed (because the dict is a defaultdict)
    """
    collapsed = {}
    for aln_id, rec in details_dict.iteritems():
        if len(rec) == 0:  # to deal with the downsides of defaultdict
            continue
        if isinstance(rec[0], list):  # hack to determine if this is a list of lists
            rec = "\n".join(["\t".join(map(str, x)) for x in rec])
        else:
            rec = "\t".join(map(str, rec))
        collapsed[aln_id] = "".join([rec, "\n"])
    return collapsed


def run_transmap_eval(cur, genome, biotype, trans_map_eval, good=True):
    """
    Convenience wrapper for getting all Good/Pass transcripts for transMap
    """
    query = trans_map_eval(genome, biotype, good)
    return get_query_ids(cur, query)


def get_stats(cur, genome, mode="transMap"):
    """
    Returns a dictionary mapping each aln_id to [aln_id, %ID, %COV].
    We replace all NULL values with zero to make things easier in consensus finding.
    """
    if mode == "transMap":
        query = ("SELECT AlignmentId,IFNULL(AlignmentCoverage, 0),IFNULL(AlignmentIdentity, 0) "
                 "FROM attributes.'{}'").format(genome)
    else:
        query = ("SELECT AugustusAlignmentId,IFNULL(AlignmentCoverage, 0),IFNULL(AlignmentIdentity, 0) "
                 "FROM augustus_attributes.'{}'").format(genome)
    return get_query_dict(cur, query)


def highest_cov_aln(cur, genome):
    """
    Returns the set of alignment IDs that represent the best alignment for each source transcript (that mapped over)
    Best is defined as highest %COV. Also reports the associated coverage and identity values.
    """
    tm_stats = get_stats(cur, genome, mode="transMap")
    combined_covs = defaultdict(list)
    for aln_id, (cov, ident) in tm_stats.iteritems():
        tx_id = psl_lib.strip_alignment_numbers(aln_id)
        combined_covs[tx_id].append([aln_id, cov, ident])
    best_cov = {}
    for tx_id, vals in combined_covs.iteritems():
        best_cov[tx_id] = sorted(vals, key=lambda x: [-x[1], -x[2]])[0]
    return best_cov


def get_highest_cov_alns(cur, genomes):
    """
    Dictionary mapping each genome to a dictionary reporting each highest coverage alignment and its metrics
    """
    return {genome: highest_cov_aln(cur, genome) for genome in genomes}
