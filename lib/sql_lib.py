"""
Convenience library for interfacting with a sqlite database.
"""
import os
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


def attach_databases(comp_ann_path, has_augustus=False):
    """
    Attaches all of the databases. Expects comp_ann_path to be the path that comparativeAnnotator wrote to.
    If has_augustus is True, expects this folder to have a augustus database.
    """
    classify_path = os.path.join(comp_ann_path, "classify.db")
    attr_path = os.path.join(comp_ann_path, "attributes.db")
    details_path = os.path.join(comp_ann_path, "details.db")
    assert all([os.path.exists(x) for x in [classify_path, attr_path, details_path]])
    con = sql.connect(classify_path)
    cur = con.cursor()
    attach_database(con, attr_path, "attributes")
    attach_database(con, details_path, "details")
    if has_augustus:
        aug_classify_path = os.path.join(comp_ann_path, "augustus_classify.db")
        aug_details_path = os.path.join(comp_ann_path, "augustus_details.db")
        aug_attributes_path = os.path.join(comp_ann_path, "augustus_attributes.db")
        assert all([os.path.exists(x) for x in [aug_classify_path, aug_details_path, aug_attributes_path]])
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


def get_biotype_ids(cur, genome, biotype, category="Transcript"):
    """
    Returns a set of aln_ids which are members of the biotype. Category should be Transcript or Gene
    """
    assert category in ["Transcript", "Gene"]
    query = "SELECT AlignmentId FROM attributes.'{}' WHERE {}Type='{}'".format(genome, category, biotype)
    return {x[0] for x in cur.execute(query).fetchall()}


def get_stats(cur, genome, category="transMap"):
    """
    Returns a dictionary mapping each aln_id to [aln_id, %ID, %COV]
    """
    if category == "transMap":
        cmd = "SELECT AlignmentId,AlignmentIdentity,AlignmentCoverage FROM attributes.'{}'".format(genome)
    else:
        cmd = "SELECT AlignmentId,AlignmentIdentity,AlignmentCoverage FROM augustus_attributes.'{}'".format(genome)
    return {x[0]: x for x in cur.execute(cmd).fetchall()}


def write_dict(data_dict, database_path, table):
    """
    Writes a dict of dicts to a sqlite database.
    """
    df = pd.DataFrame.from_dict(data_dict)
    df = df.sort_index()
    with ExclusiveSqlConnection(database_path) as con:
        df.to_sql(table, con, if_exists="replace", index_label="AlignmentId")


def collapse_details_dict(details_dict):
    """
    Collapses a details dict into a string so that the database record is a string. Properly accounts for the two
    possibilities: the record is a single BED entry (a string) or the record is a list of BED entries
    """
    collapsed = {}
    for aln_id, rec in details_dict.iteritems():
        if isinstance(rec[0], list):  # hack to determine if this is a list of lists
            rec = "\n".join(["\t".join(map(str, x)) for x in rec])
        else:
            rec = "\t".join(map(str, rec))
        collapsed[aln_id] = "".join([rec, "\n"])
    return collapsed