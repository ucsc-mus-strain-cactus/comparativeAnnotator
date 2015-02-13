"""
Convenience library for interfacting with a sqlite database.

Designed to handle concurrency issues when writing tons of stuff as fast as possible.

Author: Ian Fiddes
"""

import sqlite3 as sql
from itertools import izip, izip_longest


class ExclusiveSqlConnection(object):
    """meant to be used with a with statement to ensure proper closure"""

    def __init__(self, path, timeout=600):
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

    created without an index for speed of insertions.
    """
    cur.execute("""CREATE TABLE '{}' ({} TEXT PRIMARY KEY) WITHOUT ROWID""".format(table, primary_key))
    for n, t in columns:
        cur.execute("""ALTER TABLE '{}' ADD COLUMN {} {} """.format(table, n, t))


def numberOfRows(cur, table):
    """
    Returns the number of rows in the provided table
    """
    cur.execute("SELECT Count(*) FROM {}".format(table))
    return cur.fetchone()[0]


def insertRow(cur, table, primary_key_column, primary_key, columns, values):
    """
    Inserts a new row into the table.
    table = table name to be inserted into.
    primary_key_column = column name of primary key
    primary_key = unique row id for this row
    columns = list of columns to be inserted
    values = list of values to be inserted into columns
    """
    columns = ", ".join([primary_key_column] + columns)
    values = [primary_key] + values
    cmd = """INSERT INTO '{}' ({}) VALUES ({})""".format(table, column_string, ", ".join(["?"] * len(values)))
    cur.execute(cmd, values)


def insertRows(cur, table, primary_key_column, columns, valueIter):
    """
    Inserts n new rows in a table with executemany. Faster than one at a time.
    table = table name to be inserted into.
    columns = list of columns to be inserted
    valueIter = list of lists of values to be inserted.

    The first item of each sublist in valueIter should be the primary key for this row.
    """
    num_columns = len(columns) + 1
    column_string = ", ".join([primary_key_column] + columns)
    cmd = """INSERT INTO '{}' ({}) VALUES ({})""".format(table, column_string, ", ".join(["?"] * num_columns))
    cur.executemany(cmd, valueIter)


def updateRow(cur, table, primary_key_column, primary_key, col_to_change, value):
    """
    Updates a row keyed by primary_key to change the value at col_to_change to
    the passed value.
    """
    cmd = """UPDATE '{}' SET {}=? WHERE {}=?""".format(table, col_to_change,
            primary_key_column)
    cur.execute(cmd, (value, primary_key))


def updateRows(cur, table, primary_key_column, col_to_change, valueIter):
    """
    Updates a row keyed by primary_key to change the value at col_to_change to
    the passed value.
    """
    cmd = """UPDATE '{}' SET {}=? WHERE {}=?""".format(table, col_to_change,
            primary_key_column)
    cur.executemany(cmd, valueIter)


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
    updateRow(cur, table, primary_key_column, primary_key, col_to_change, value)


def attachDatabase(con, path, name):
    """
    Attaches another database found at path to the name given in the given connection.
    """
    con.execute("ATTACH DATABASE '{}' AS {}".format(path, name))


def selectBetweenDatabasesWithAnd(cur, altName, returnColumn, columns, values, primaryKey, table):
    """
    Selects a column from the alternate database defined by altName that fit criteria in
    the primary database (main). columns should be a list of columns we are selecting on
    and values a matching list of values to check for. Returns all entries in altName that
    match this select statement. returnColumn is the column whose values in altName we want.

    Assumes that attachDatabase() as been run on a connection so that altName is an attached database.
    table should exist in both databases and is the table we are joining.
    """
    cmd = ("""SELECT {altName}.{table}.{returnColumn} FROM main.{table} JOIN {altName}.{table} USING ('{primaryKey}')"""
            """ WHERE {altName}.{table}.{returnColumn} IS NOT NULL AND main.{table}.{primaryKey} = """
            """{altName}.{table}.{primaryKey}""").format(altName=altName, returnColumn=returnColumn, table=table,
                                                         primaryKey=primaryKey)
    for col in columns:
        cmd += " AND main.{}.{} = ?".format(table, col)
    print cmd
    q = cur.execute(cmd, values)
    return q.fetchall()


def selectBetweenDatabasesWithOr(cur, altName, returnColumn, columns, values, primaryKey, table):
    """
    Selects a column from the alternate database defined by altName that fit criteria in
    the primary database (main). columns should be a list of columns we are selecting on
    and values a matching list of values to check for. Returns all entries in altName that
    match this select statement. returnColumn is the column whose values in altName we want.

    Assumes that attachDatabase() as been run on a connection so that altName is an attached database.
    table should exist in both databases and is the table we are joining.
    """
    cmd = ("""SELECT {altName}.{table}.{returnColumn} FROM main.{table} JOIN {altName}.{table} USING ('{primaryKey}')"""
            """ WHERE {altName}.{table}.{returnColumn} IS NOT NULL AND main.{table}.{primaryKey} = """
            """{altName}.{table}.{primaryKey}""").format(altName=altName, returnColumn=returnColumn, table=table,
                                                         primaryKey=primaryKey)
    for col in columns:
        cmd += " OR main.{}.{} = ?".format(table, col)
    print cmd
    q = cur.execute(cmd, values)
    return q.fetchall()


def selectBetweenDatabases(cur, altName, returnColumn, columns, values, modifiers, primaryKey, table):
    """
    Selects a column from the alternate database defined by altName that fit criteria in
    the primary database (main). columns should be a list of columns we are selecting on
    and values a matching list of values to check for. Returns all entries in altName that
    match this select statement. returnColumn is the column whose values in altName we want.

    modifiers should say either "AND" or "OR" describing how each comparison in the sequence should be join

    Assumes that attachDatabase() as been run on a connection so that altName is an attached database.
    table should exist in both databases and is the table we are joining.
    """
    cmd = ("""SELECT {altName}.{table}.{returnColumn} FROM main.{table} JOIN {altName}.{table} USING ('{primaryKey}')"""
            """ WHERE {altName}.{table}.{returnColumn} IS NOT NULL AND main.{table}.{primaryKey} = """
            """{altName}.{table}.{primaryKey}""").format(altName=altName, returnColumn=returnColumn, table=table,
                                                         primaryKey=primaryKey)
    for col, mod in izip(columns, modifiers):
        cmd += " {} main.{}.{} = ?".format(mod, table, col)
    print cmd
    q = cur.execute(cmd, values)
    return q.fetchall()
