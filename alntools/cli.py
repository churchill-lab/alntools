# -*- coding: utf-8 -*-
from __future__ import print_function

import glob
import os

import click

from . import bin_utils
from . import methods
from . import utils
from . import viewer
from . import __logo_text__, __version__


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, message=__logo_text__)
def cli():
    """
    alntools

    Simple tools for alignment

    """


@cli.command('split', options_metavar='<options>', short_help='split a BAM file into smaller bam files')
@click.argument('bam_file', metavar='bam_file', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.argument('number', metavar='number', type=int)
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=True), help="output directory")
@click.option('-v', '--verbose', count=True, help='use \'-vv\' or \'-vvv\' for the more output')
def split(bam_file, number, directory, verbose):
    """
    Split a BAM file into smaller bam files.
    """
    utils.configure_logging(verbose)
    methods.split_bam(bam_file, number, directory)


@cli.command('bam2ec', options_metavar='<options>', short_help='convert a BAM file to EC')
@click.argument('bam_file', metavar='bam_file', type=click.Path(exists=True, resolve_path=True, dir_okay=True))
@click.argument('ec_file', metavar='ec_file', type=click.Path(resolve_path=True, dir_okay=False, writable=True))
@click.option('-c', '--chunks', default=0, help="number of chunks to process")
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=True), help="temp directory")
@click.option('-m', '--mincount', default=2000, help="minimum count")
@click.option('--multisample', is_flag=True)
@click.option('-p', '--number_processes', default=-1, help="number of processes")
@click.option('--rangefile', type=click.Path(exists=False, resolve_path=True, file_okay=True, dir_okay=False, writable=True), help="range file")
@click.option('-s', '--sample', help="sample identifier")
@click.option('-t', '--targets', metavar='FILE', type=click.Path(exists=True, resolve_path=True, file_okay=True, dir_okay=False), help="target file")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def bam2ec(bam_file, ec_file, chunks, directory, mincount, multisample, number_processes, rangefile, sample, targets, verbose):
    """
    Convert a BAM file (bam_file) to a binary EC file (ec_file).
    """
    utils.configure_logging(verbose)
    if multisample:
        if sample:
            print('-s, --sample should NOT be specified with --multisample')
            return
        methods.bam2ec_multisample(bam_file, ec_file, chunks, mincount, directory, number_processes, rangefile, targets)
    else:
        methods.bam2ec(bam_file, ec_file, chunks, directory, number_processes, rangefile, sample, targets)


@cli.command('bam2emase', options_metavar='<options>', short_help='convert a BAM file to APM')
@click.argument('bam_file', metavar='bam_file', type=click.Path(exists=True, resolve_path=True, dir_okay=True))
@click.argument('emase_file', metavar='emase_file', type=click.Path(resolve_path=True, dir_okay=False, writable=True))
@click.option('-c', '--chunks', default=0, help="number of chunks to process")
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=True), help="temp directory")
@click.option('-m', '--mincount', default=2000, help="minimum count")
@click.option('--multisample', is_flag=True)
@click.option('-p', '--number_processes', default=-1, help="number of processes")
@click.option('--rangefile', type=click.Path(exists=False, resolve_path=True, file_okay=True, dir_okay=False, writable=True), help="range file")
@click.option('-t', '--targets', metavar='FILE', type=click.Path(exists=True, resolve_path=True, file_okay=True, dir_okay=False), help="target file")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def bam2emase(bam_file, emase_file, chunks, directory, mincount, multisample, number_processes, rangefile, targets, verbose):
    """
    Convert a BAM file (bam_file) to an EMASE file (emase_file).
    """
    utils.configure_logging(verbose)

    if multisample:
        methods.bam2emase_multisample(bam_file, emase_file, chunks, mincount, directory, number_processes, rangefile, targets)
    else:
        methods.bam2emase(bam_file, emase_file, chunks, directory, number_processes, rangefile, targets)


@cli.command('bam2both', options_metavar='<options>', short_help='convert a BAM file to EC and Emase')
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
def bam2both(bam_file, ec_file, emase_file, chunks, directory, mincount, multisample, number_processes, rangefile, sample, verbose, targets):
    """
    Convert a BAM file (bam_file) to a binary EC file (ec_file) and an EMASE file (emase_file).
    """
    utils.configure_logging(verbose)
    if multisample:
        if sample:
            print('-s, --sample should NOT be specified with --multisample')
            return
        methods.bam2both_multisample(bam_file, ec_file, emase_file, chunks, mincount, directory, number_processes, rangefile, targets)
    else:
        methods.bam2both(bam_file, ec_file, emase_file, chunks, directory, number_processes, rangefile, sample, targets)


@cli.command('salmon2ec', options_metavar='<options>', short_help='convert a salmon eq_classes file to EC')
@click.argument('salmon_dir', metavar='salmon_dir', type=click.Path(exists=True, resolve_path=True, dir_okay=True))
@click.argument('ec_file', metavar='ec_file', type=click.Path(resolve_path=True, dir_okay=False, writable=True))
@click.option('-s', '--sample', metavar='sample', default='NA', help="sample identifier")
@click.option('-t', '--targets', metavar='FILE', type=click.Path(exists=True, resolve_path=True, file_okay=True, dir_okay=False), help="target file")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def salmon2ec(salmon_dir, ec_file, sample, targets, verbose):
    """
    Convert a salmon eq_classes file to an EC file (ec_file).
    """
    utils.configure_logging(verbose)
    methods.salmon2ec(salmon_dir, ec_file, sample, targets)


@cli.command('range', options_metavar='<options>', short_help='test')
@click.argument('input', nargs=-1, type=click.Path(resolve_path=True, dir_okay=True))
@click.argument('range_file', metavar='range_file', type=click.Path(exists=False, resolve_path=True, file_okay=True, dir_okay=False, writable=True))
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=True), help="temp directory")
@click.option('-t', '--targets', metavar='FILE', type=click.Path(exists=True, resolve_path=True, file_okay=True, dir_okay=False), help="target file")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def range(input, range_file, targets, directory, verbose):
    """
    Create range file for BAM files
    """
    utils.configure_logging(verbose)
    files = utils.get_bam_files(input)
    methods.generate_bam_ranges(files, range_file, targets, directory)


@cli.command('emase2db_configure', options_metavar='<options>', short_help='generate sample file for viewer')
@click.argument('sample_file', metavar='sample_file', type=click.Path(resolve_path=True, dir_okay=False))
@click.option('-d', '--directory', type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True, writable=False), help="top level input directory")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def emase2db_config(sample_file, directory, verbose):
    """
    Generate config file for viewer.
    """
    utils.configure_logging(verbose)
    methods.emase2db_config(sample_file, directory)


@cli.command('emase2db', options_metavar='<options>', short_help='generate database file for viewer')
@click.argument('sample_file', metavar='sample_file', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.argument('gene_file', metavar='gene_file', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.argument('db_file', metavar='db_file', type=click.Path(resolve_path=True, dir_okay=False))
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def emase2db(sample_file, gene_file, db_file, verbose):
    """
    Generate database file for viewer
    """
    utils.configure_logging(verbose)
    methods.emase2db(sample_file, gene_file, db_file)


@cli.command('view', options_metavar='<options>', short_help='generate database file for viewer')
@click.argument('db_file', metavar='db_file', type=click.Path(resolve_path=True, dir_okay=False))
@click.option('-p', '--port', default=8888, help="port number")
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def v(db_file, port, verbose):
    """
    Generate database file for viewer
    """
    utils.configure_logging(verbose)
    viewer.start(db_file, port)


@cli.command('dumpec', options_metavar='<options>', short_help='dumpe file information')
@click.argument('ec_file', metavar='ecfile', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def dumpec(ec_file, verbose):
    """
    Dump information about an EC file (ec_file).
    """
    utils.configure_logging(verbose)
    methods.dumpec(ec_file)


@cli.command('merge', options_metavar='<options>', short_help='merge multiple ec files')
@click.option('-i', '--input',
              metavar='input',
              type=click.Path(exists=True, resolve_path=True, dir_okay=False),
              multiple=True,
              help="input file, can specify multiple")
@click.option('-d', '--directory',
              metavar='directory',
              type=click.Path(exists=True, resolve_path=True, file_okay=False, dir_okay=True),
              help="input directory")
@click.option('-o', '--output',
              metavar='output',
              type=click.Path(resolve_path=True, dir_okay=False),
              help="output file")
@click.option('-v', '--verbose',
              count=True,
              help='the more times listed, the more output')
def merge(input, directory, output, verbose):
    """
    Combine multiple ec files.
    """
    utils.configure_logging(verbose)
    #methods.dumpec(ec_file)
    #print(input)
    #print(directory)
    #print(output)

    input_files = list(input)

    if directory:
        bin_files = glob.glob(os.path.join(directory, "*.*"))
        if len(bin_files) == 0:
            print('No bin files found in directory: {}'.format(directory))
            return None

        if input_files is None:
            input_files = bin_files
        else:
            input_files.extend(bin_files)


    #print(input_files)
    bin_utils.combine(input_files, output)


@cli.command('ec2emase', options_metavar='<options>', short_help='convert ec file to APM')
@click.argument('ec_file', metavar='ec_file', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.argument('apm_file', metavar='apm_file', type=click.Path(resolve_path=True, dir_okay=False))
@click.option('-v', '--verbose', count=True, help='the more times listed, the more output')
def emase2db(ec_file, apm_file, verbose):
    """
    Generate database file for viewer
    """
    utils.configure_logging(verbose)
    methods.ec2apm(ec_file, apm_file)


if __name__ == '__main__':
    cli()

