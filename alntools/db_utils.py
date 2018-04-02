# -*- coding: utf-8 -*-
from collections import OrderedDict

import sqlite3

expected_gene_info_header = ['GeneID', 'Symbol', 'Chr', 'Start', 'End', 'Strand']
expected_gene_edges_header = ['GeneID1', 'GeneID2', 'Weight']


def _dictify_row(cursor, row):
    """Turns the given row into a dictionary where the keys are the column names"""
    d = OrderedDict()
    for i, col in enumerate(cursor.description):
        d[col[0]] = row[i]
    return d


def dictify_cursor(cursor):
    """converts all cursor rows into dictionaries where the keys are the column names"""
    return (_dictify_row(cursor, row) for row in cursor)


def init_db(db_file):
    db_con = sqlite3.connect(db_file)

    try:
        c = db_con.cursor()

        c.execute(
            '''
            CREATE TABLE gene_count_totals (
                sample_id   TEXT NOT NULL,
                gene_id     TEXT NOT NULL,
                read_count  REAL NOT NULL,
                PRIMARY KEY(sample_id, gene_id))
            '''
        )

        c.execute(
            '''
            CREATE TABLE gene_edge (
                sample_id   TEXT  NOT NULL,
                gene_id_1   TEXT  NOT NULL,
                gene_id_2   TEXT  NOT NULL,
                read_count  REAL  NOT NULL,
                PRIMARY KEY(sample_id, gene_id_1, gene_id_2))
            '''
        )

        c.execute(
            '''
            CREATE TABLE gene_info (
                gene_id       TEXT PRIMARY KEY  NOT NULL,
                chromosome    TEXT              NOT NULL,
                start_pos_bp  INTEGER           NOT NULL,
                end_pos_bp    INTEGER           NOT NULL,
                strand        TEXT              NOT NULL,
                symbol        TEXT              NOT NULL,
                name          TEXT              NULL)
            '''
        )
        c.execute('''CREATE INDEX gene_info_symbol_idx ON gene_info (symbol)''')

        c.execute(
            '''
            CREATE TABLE sample (
                sample_id TEXT PRIMARY KEY NOT NULL,
                sample_file_name TEXT NOT NULL)
            ''')

        db_con.commit()
    except:
        db_con.rollback()
        raise

    return db_con


def non_zero_expr_intervals(db_con, sample_id):
    # TODO not using sample_id for the moment
    c = db_con.cursor()
    c.execute(
        '''
        SELECT chromosome, start_pos_bp, end_pos_bp, read_count
        FROM gene_info INNER JOIN gene_count_totals ON gene_info.gene_id = gene_count_totals.gene_id
        WHERE read_count>0.0 AND sample_id=?
        ''',
        (sample_id, ))
    return list(dictify_cursor(c))


def multiread_counts(db_con, sample_id, gene_id):
    c = db_con.cursor()
    c.execute(
        '''
        SELECT gene_info.*, gene_edge.read_count
        FROM gene_edge INNER JOIN gene_info ON gene_edge.gene_id_2 = gene_info.gene_id
        WHERE sample_id=? AND gene_id_1=?
        ORDER BY gene_edge.read_count DESC
        ''',
        (sample_id, gene_id))
    return list(dictify_cursor(c))


def multiread_counts_all_genes(db_con, sample_id):
    c = db_con.cursor()
    c.execute(
        '''
        SELECT g1_info.*, g2_info.*, gene_edge.read_count
        FROM
          gene_edge
          INNER JOIN gene_info AS g1_info ON gene_edge.gene_id_1 = g1_info.gene_id
          INNER JOIN gene_info AS g2_info ON gene_edge.gene_id_2 = g2_info.gene_id
        WHERE sample_id=?
        ORDER BY g1_info.gene_id, g2_info.gene_id
        ''',
        (sample_id, ))
    return list(dictify_cursor(c))


def gene_info(db_con, gene_id):
    c = db_con.cursor()
    c.execute('''SELECT * FROM gene_info WHERE gene_id=?''', (gene_id, ))
    rows = list(dictify_cursor(c))
    return rows[0] if len(rows) == 1 else None


def all_gene_info(db_con):
    c = db_con.cursor()
    c.execute('''SELECT * FROM gene_info''')
    return list(dictify_cursor(c))


def lookup_gene_id(db_con, gene_id):
    c = db_con.cursor()
    c.execute('''SELECT gene_id FROM gene_info WHERE UPPER(symbol) = UPPER(?)''', (gene_id, ))
    ids = [row[0] for row in c]
    if len(ids) >= 1:
        # TODO this isn't really the right thing to do in the case of >1
        return ids[0]
    else:
        return None


def all_sample_ids(db_con):
    c = db_con.cursor()
    c.execute('''SELECT * from sample''')
    return list(dictify_cursor(c))



def add_gene_count_total(cursor, sample_id, gene_id, read_count):
    cursor.execute('''INSERT INTO gene_count_totals VALUES (?, ?, ?)''', (sample_id, gene_id, read_count))


def add_gene_edge(cursor, sample_id, gene_id_1, gene_id_2, read_count):
    cursor.execute('''INSERT INTO gene_edge VALUES (?, ?, ?, ?)''', (sample_id, gene_id_1, gene_id_2, read_count))


def add_gene_info(cursor, gene_id, chromosome, start_pos_bp, end_pos_bp, strand, symbol, name):
    cursor.execute(
        '''INSERT INTO gene_info VALUES (?, ?, ?, ?, ?, ?, ?)''',
        (gene_id, chromosome, start_pos_bp, end_pos_bp, strand, symbol, name))


def add_sample_id(cursor, sample_id, sample_file_name):
    cursor.execute(
        '''INSERT INTO sample VALUES (?, ?)''', (sample_id, sample_file_name))



