"""
Convenience library for interfacting with a sqlite database. Designed to handle concurrency
issues when writing tons of stuff.

Author: Ian Fiddes
"""

import sqlite3 as sql


class ExclusiveSqlConnection(object):
    """meant to be used with a with statement to ensure proper closure"""

    def __init__(self, path, timeout=6000):
        self.path = path
        self.timeout = timeout

    def __enter__(self):
        self.con = sql.connect(self.path, timeout = self.timeout, isolation_level = "EXCLUSIVE")
        try:
            self.con.execute("BEGIN EXCLUSIVE")
        except sql.OperationalError:
            print ("Database still locked after {} seconds.".format(self.timeout))
        return self.con.cursor()

    def __exit__(self, exception_type, exception_val, trace):
        self.con.commit()
        self.con.close()


def hasTable(cur, table):
    """checks to make sure this sql database has a specific table"""
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='table_name'")
    rows = cur.fetchall()
    if table in rows:
        return True
    else:
        return False


def updateTable(cur, table, set_col, set_val, where_col, where_val):
    """
    table = name of table to be updated
    set_col = column where change should occur
    set_val = value in that column to change
    where_col = column whose value determines if changes occur
    where_val = value in where_col that determines which rows get changed

    This breaks a big no-no of SQL and allows injection attacks. But, who cares?
    This isn't a web application.
    """
    cur.execute("UPDATE {} SET {}={} WHERE {}={}".format(table, set_col, set_val,
            where_col, where_val))


def initializeTable(cur, table, columns, primary_key):
    """
    Builds a empty <table> in the sqlite db opened by <cur> with <[columns]>
    columns should be a list of name, type pairs
    primary_key should be the column name that will be the primary key.

    This breaks a big no-no of SQL and allows injection attacks. But, who cares?
    This isn't a web application.
    """
    cur.execute("""CREATE TABLE '{}' ({} TEXT PRIMARY KEY)""".format(table, primary_key))
    for n, t in columns:
        cur.execute("""ALTER TABLE '{}' ADD COLUMN {} {} """.format(table, n, t))


def numberOfRows(cur, table):
    """
    Returns the number of rows in the provided table
    """
    cur.execute("SELECT Count(*) FROM {}".format(table))
    return cur.fetchone()[0]


def insertRow(cur, table, primary_key_column, primary_key):
    cmd = """INSERT INTO '{}' ({}) VALUES ('{}')""".format(table, primary_key_column, primary_key)
    cur.execute(cmd)


def upsert(cur, table, primary_key_column, primary_key, col_to_change, value):
    """
    Wrapper for a 'upsert' command (which sqlite lacks). Given a table and a primary key
    will insert a new row with this primary key as necessary. Current implementation
    only inserts a value for one column.

    This breaks a big no-no of SQL and allows injection attacks. But, who cares?
    This isn't a web application.
    """
    #TODO: make col_to_change and value a list and iterate over them?
    cmd = """INSERT OR IGNORE INTO '{}' ({}, {}) VALUES (?, ?)""".format(table, 
            primary_key_column, col_to_change)
    cur.execute(cmd, (primary_key, value))
    cmd = """UPDATE '{}' SET {}=? WHERE {}=?""".format(table, col_to_change,
            primary_key_column)
    cur.execute(cmd, (value, primary_key))
