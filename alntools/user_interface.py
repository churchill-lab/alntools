# -*- coding: utf-8 -*-
from six import iteritems

import flask
import math
import os
import sqlite3

from flask import request

from alntools import __logo_text__, __version__
from alntools.utils import str2bool, configure_logging

source_dir = os.path.dirname(os.path.realpath(__file__))
app = flask.Flask(__name__, static_folder=os.path.join(source_dir, 'static'), static_url_path='')





@app.route('/')
def index():
    return flask.render_template('index.html', logo_text=__logo_text__, version=__version__)



@app.route('/jqueryfiletree', methods=['GET', 'POST'])
def jqueryfiletree():
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")
    #print(request.values)

    only_folders = str2bool(request.values.get("onlyFolders", "0"))
    only_files = str2bool(request.values.get("onlyFiles", "0"))
    checkbox = str2bool(request.values.get("multiSelect", ""))

    if checkbox:
        checkbox = "<input value='{}' isdir='{}' type='checkbox'/>"
    else:
        checkbox = ""

    r = ['<ul class="jqueryFileTree" style="display: none;">']
    try:
        r = ['<ul class="jqueryFileTree" style="display: none;">']
        d = request.values.get("dir", "/")

        for f in os.listdir(d):
            ff = os.path.join(d, f)
            if os.path.isdir(ff) and (not only_files or only_folders):
                checkbox_formatted = checkbox.format(ff, 'true')
                r.append(
                    '<li class="directory collapsed">' + checkbox_formatted + ' <a rel="%s/">%s</a></li>' % (
                    ff, f))
            elif (not only_folders) or only_files:
                e = os.path.splitext(f)[1]  # get .ext and remove dot
                if len(e) > 0:
                    e = e[1:]
                checkbox_formatted = checkbox.format(ff, 'false')
                r.append('<li class="file ext_{}">{} <a rel="{}">{}</a></li>'.format(e, checkbox_formatted, ff, f))
        r.append('</ul>')
    except Exception, e:
        r.append('Could not load directory: %s' % str(e))
    r.append('</ul>')

    return ''.join(r)


@app.route('/bam2both')
def bam2both_web():
    return flask.render_template('bam2both.html', logo_text=__logo_text__, version=__version__)


@app.route('/bam2both/submit')
def bam2both_submit():
    only_folders = str2bool(request.values.get("onlyFolders", "0"))
    only_files = str2bool(request.values.get("onlyFiles", "0"))
    checkbox = str2bool(request.values.get("multiSelect", ""))

    configure_logging(verbose)
    if multisample:
        if sample:
            print('-s, --sample should NOT be specified with --multisample')
            return
        methods.bam2both_multisample(bam_file, ec_file, emase_file, chunks, mincount, directory, number_processes, rangefile, targets)
    else:
        methods.bam2both(bam_file, ec_file, emase_file, chunks, directory, number_processes, rangefile, sample, targets)

'''

    @click.argument('bam_file', metavar='bam_file', type=click.Path(exists=True, resolve_path=True, dir_okay=True))
@click.argument('ec_file', metavar='ec_file', type=click.Path(resolve_path=True, dir_okay=False, writable=True))
@click.argument('emase_file', metavar='emase_file', type=click.Path(resolve_path=True, dir_okay=False, writable=True))
@click.option('-c', '--chunks', default=0, help="number of chunks to process")
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=True), help="temp directory")
@click.option('-m', '--mincount', default=2000, help="minimum count")
@click.option('--multisample', is_flag=True)
@click.option('-p', '--number_processes', default=-1, help="number of processes")
@click.option('--rangefile', type=click.Path(exists=False, resolve_path=True, file_okay=True, dir_okay=False, writable=True), help="range file")
@click.option('-s', '--sample', help="sample identifier")
@click.option('-t', '--targets', metavar='FILE', type=click.Path(exists=True, resolve_path=True, file_okay=True, dir_okay=False), help="target file")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
'''



def start(port):
    try:
        import webbrowser
        webbrowser.open("http://localhost:{}/".format(port))
    except:
        pass

    app.debug = True
    app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False
    app.run(host='0.0.0.0', port=port)


