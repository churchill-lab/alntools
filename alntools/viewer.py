import flask
import math
import db_utils
import os
import sqlite3

source_dir = os.path.dirname(os.path.realpath(__file__))
app = flask.Flask(__name__, static_folder=os.path.join(source_dir, 'static'), static_url_path='')

DB_FILE = None


def make_db_connection():
    return sqlite3.connect(DB_FILE)


# example:
# http://roo:8888/log2expr-for-sample-sample_id/bin-size-1000000bp.json
@app.route('/log2expr-for-sample-<sample_id>/bin-size-<int:bin_size_bp>bp.json')
def binned_log2_expr_json(sample_id, bin_size_bp):
    expr_intervals = db_utils.non_zero_expr_intervals(make_db_connection(), sample_id)
    max_expr = 0
    expr_map = dict()
    bin_count_per_chr = dict()

    for expr_interval in expr_intervals:
        start_bin_index = int(expr_interval['start_pos_bp']) / bin_size_bp
        end_bin_index = int(expr_interval['end_pos_bp']) / bin_size_bp
        chr = expr_interval['chromosome']
        log2_expr = math.log(expr_interval['read_count'], 2)

        if log2_expr > max_expr:
            max_expr = log2_expr

        for bin_index in xrange(start_bin_index, end_bin_index + 1):
            if chr not in bin_count_per_chr or bin_index + 1 > bin_count_per_chr[chr]:
                bin_count_per_chr[chr] = bin_index + 1

            if chr in expr_map:
                chr_bin_dict = expr_map[chr]
            else:
                chr_bin_dict = dict()
                expr_map[chr] = chr_bin_dict

            if bin_index not in chr_bin_dict or log2_expr > chr_bin_dict[bin_index]:
                chr_bin_dict[bin_index] = log2_expr

    chr_expr_bins = dict()
    for chr, bin_count in bin_count_per_chr.iteritems():
        log2_expr_bins = [0] * bin_count
        for bin_index, bin_expr in expr_map[chr].iteritems():
            log2_expr_bins[bin_index] = bin_expr
        chr_expr_bins[chr] = log2_expr_bins

    return flask.jsonify(
        max_expr=max_expr,
        chr_expr_bins=chr_expr_bins,
        bin_size_bp=bin_size_bp,
    )


# example:
# http://roo:8888/multiread-weights-HARD_CODED_SAMPLE_ID/gene-ENSMUSG00000098192.json
@app.route('/multiread-weights-<sample_id>/gene-<gene_id>.json')
def multiread_weights_json(sample_id, gene_id):
    db_con = make_db_connection()

    ens_gene_id = gene_id
    query_gene_info = db_utils.gene_info(db_con, ens_gene_id)
    if query_gene_info is None:
        ens_gene_id = db_utils.lookup_gene_id(db_con, gene_id)

    if ens_gene_id is None:
        raise 'Failed to find gene named: {}'.format(gene_id)

    weights = db_utils.multiread_counts(db_con, sample_id, ens_gene_id)
    query_gene_info = db_utils.gene_info(db_con, ens_gene_id)
    return flask.jsonify(
        sample_id=sample_id,
        query_gene_info=query_gene_info,
        multiread_genes=weights,
    )


# @app.route('/multiread-weights-<sample_id>/all-genes.json')
# def all_multiread_weights_json(sample_id):
#     db_con = make_db_connection()
#
#     ens_gene_id = gene_id
#     query_gene_info = multireaddb.gene_info(db_con, ens_gene_id)
#     if query_gene_info is None:
#         ens_gene_id = multireaddb.lookup_gene_id(db_con, gene_id)
#
#     if ens_gene_id is None:
#         raise 'Failed to find gene named: {}'.format(gene_id)
#
#     weights = multireaddb.multiread_counts(db_con, sample_id, ens_gene_id)
#     query_gene_info = multireaddb.gene_info(db_con, ens_gene_id)
#     return flask.jsonify(
#         sample_id=sample_id,
#         query_gene_info=query_gene_info,
#         multiread_genes=weights,
#     )


@app.route('/all-gene-info.json')
def all_gene_info():
    db_con = make_db_connection()
    gene_info = db_utils.all_gene_info(db_con)
    return flask.jsonify(gene_info=gene_info)


@app.route('/all-samples.json')
def all_samples():
    db_con = make_db_connection()
    return flask.jsonify(sample_ids=db_utils.all_sample_ids(db_con))


def start(db_file, port):
    global DB_FILE
    DB_FILE = db_file

    import webbrowser
    webbrowser.open("http://localhost:{}/multireadview.html".format(port))

    app.debug = True
    app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False
    app.run(host='0.0.0.0', port=port)

