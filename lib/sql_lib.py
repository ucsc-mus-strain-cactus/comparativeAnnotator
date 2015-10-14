"""
Convenience library for interfacting with a sqlite database.
"""
import os
import re
import sqlite3 as sql
import pandas as pd

__author__ = "Ian Fiddes"


class ExclusiveSqlConnection(object):
    """meant to be used with a with statement to ensure proper closure"""

    def __init__(self, path, timeout=600):
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


def attach_databases(comp_ann_path, mode):
    """
    Attaches all of the databases needed for this execution mode.
    """
    assert mode in ["reference", "augustus", "transMap"]
    classify_path = os.path.join(comp_ann_path, "classify.db")
    details_path = os.path.join(comp_ann_path, "details.db")
    con = sql.connect(classify_path)
    cur = con.cursor()
    attach_database(con, details_path, "details")
    if mode in ["transMap", "augustus"]:
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


def get_query_ids(cur, query):
    """
    Returns a set of aln_ids which are OK based on the definition of OK in config.py that made this query.
    In other words, expects a query of the form SELECT AlignmentId FROM stuff
    """
    return {x[0] for x in cur.execute(query).fetchall()}


def get_query_dict(cur, query):
    """
    Returns a set of aln_ids which are OK based on the definition of OK in config.py that made this query.
    In other words, expects a query of the form SELECT AlignmentId,<other stuff> FROM stuff
    """
    return {x[0]: x[1:] for x in cur.execute(query).fetchall()}


def get_biotype_aln_ids(cur, genome, biotype, mode="Transcript"):
    """
    Returns a set of aln_ids which are members of the biotype. Category should be Transcript or Gene
    """
    assert mode in ["Transcript", "Gene"]
    query = "SELECT AlignmentId FROM attributes.'{}' WHERE {}Type='{}'".format(genome, mode, biotype)
    return get_query_ids(cur, query)


def get_biotype_ids(cur, genome, biotype, mode="Transcript"):
    """
    Returns a set of aln_ids which are members of the biotype. Category should be Transcript or Gene
    """
    assert mode in ["Transcript", "Gene"]
    query = "SELECT TranscriptId FROM attributes.'{}' WHERE {}Type='{}'".format(genome, mode, biotype)
    return get_query_ids(cur, query)


def filter_biotype_ids(cur, genome, biotype, filter_set, mode="Transcript"):
    """
    Runs get_biotype_ids() filtering the resulting set based on the blacklist filter_set
    """
    ids = get_biotype_ids(cur, genome, biotype, mode)
    return ids - filter_set


def get_stats(cur, genome, mode="transMap"):
    """
    Returns a dictionary mapping each aln_id to [aln_id, %ID, %COV]
    """
    if mode == "transMap":
        query = "SELECT AlignmentId,AlignmentIdentity,AlignmentCoverage FROM attributes.'{}'".format(genome)
    else:
        query = "SELECT AlignmentId,AlignmentIdentity,AlignmentCoverage FROM augustus_attributes.'{}'".format(genome)
    return get_query_dict(cur, query)


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


def get_all_biotypes(cur, ref_genome, gene_level=False):
    """
    Returns all biotypes in the attribute database.
    """
    r = "Gene" if gene_level else "Transcript"
    query = "SELECT {}Type FROM attributes.'{}'".format(r, ref_genome)
    return get_query_ids(cur, query)


def run_transmap_eval(cur, genome, biotype, trans_map_eval, good=True):
    """
    Convenience wrapper for getting all Good/Pass transcripts for transMap
    """
    query = trans_map_eval(genome, biotype, good)
    return get_query_ids(cur, query)


def highest_cov_aln(cur, genome):
    """
    Returns the set of alignment IDs that represent the best alignment for each source transcript (that mapped over)
    Best is defined as highest %COV. Also reports the associated coverage value.
    """
    tm_stats = get_stats(cur, genome, category="transMap")
    combined_covs = defaultdict(list)
    for aln_id, (ident, cov) in tm_stats.iteritems():
        tx_id = psl_lib.strip_alignment_numbers(aln_id)
        combined_covs[tx_id].append([aln_id, ident, cov])
    best_cov = {}
    for tx_id, vals in combined_covs.iteritems():
        best_cov[tx_id] = sorted(vals, key=lambda x: -x[2])[0]
    return best_cov


def get_ids_by_chromosome(cur, genome, chromosomes=["Y", "chrY"]):
    """
    Returns all transcript IDs in the attribute database whose source chromosome is in chromomsomes
    """
    base_query = "SELECT AlignmentId FROM attributes.'{}' WHERE {}"
    filt = ["SourceChrom = '{}'".format(x) for x in chromosomes]
    filt = " OR ".join(filt)
    query = base_query.format(genome, filt)
    return get_query_ids(cur, query)