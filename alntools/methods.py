# -*- coding: utf-8 -*-
from collections import OrderedDict
from struct import pack

import csv
import gzip
import os
import re
import tempfile
import time

import numpy as np

from scipy.sparse import coo_matrix, diags

from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from . import bam_utils
from . import bam_utils_multisample
from . import barcode_utils
from . import bin_file
from . import bin_utils
from . import db_utils
from . import utils

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


def bam2ec(bam_filename, ec_filename, chunks=0, directory=None, number_processes=-1, range_filename=None, sample=None):
    bam_utils.convert(bam_filename, ec_filename, None, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename, sample=sample)


def bam2emase(bam_filename, emase_filename, chunks=0, directory=None, number_processes=-1, range_filename=None):
    bam_utils.convert(bam_filename, None, emase_filename, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2both(bam_filename, ec_filename, emase_filename, chunks=0, directory=None, number_processes=-1, range_filename=None, sample=None):
    bam_utils.convert(bam_filename, ec_filename, emase_filename, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename, sample=sample)


def bam2ec_multisample(bam_filename, ec_filename, chunks=0, minimum_count=-1, directory=None, number_processes=-1, range_filename=None):
    bam_utils_multisample.convert(bam_filename, ec_filename, None, num_chunks=chunks, minimum_count=minimum_count, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2both_multisample(bam_filename, ec_filename, emase_filename, chunks=0, minimum_count=-1, directory=None, number_processes=-1, range_filename=None):
    bam_utils_multisample.convert(bam_filename, ec_filename, emase_filename, num_chunks=chunks, minimum_count=minimum_count, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def split_bam(bam_filename, num_chunks, directory=None):
    bam_utils.split_bam(bam_filename, num_chunks, directory)


def generate_bam_ranges(input_files, range_filename, target_filename=None, temp_dir=None):
    bam_utils.generate_bam_ranges(input_files, range_filename, target_filename, temp_dir)


def emase2db_config(sample_file, directory=None):
    start_time = time.time()

    working_directory = os.getcwd() if directory is None else directory
    working_directory = os.path.abspath(working_directory)

    LOG.info('Generating sample file: {}'.format(sample_file))
    LOG.info('Top level directory: {}'.format(working_directory))

    samples = {}
    all_files = {}

    try:
        for root, dirs, files in os.walk(working_directory):
            for file in files:
                if file.endswith(".h5"):
                    found_file = os.path.join(root, file)
                    LOG.debug(found_file)

                    sample_name = file

                    if sample_name in samples:
                        samples[sample_name] += 1
                        sample_name = "{}_{}".format(sample_name, samples[sample_name])
                    else:
                        samples[sample_name] = 1

                    all_files[sample_name] = found_file

        with open(sample_file, "w") as out_fd:
            for k, v in all_files.iteritems():
                LOG.debug("{}\t{}".format(k, v))
                out_fd.write("{}\t{}\n".format(k, v))

    except Exception as e:
        LOG.error('Error: {}'.format(str(e)))

    LOG.info("Config file created in {}".format(utils.format_time(start_time, time.time())))


def get_num_shared_multireads(alnmat):
    hapsum = alnmat.sum(axis=APM.Axis.HAPLOTYPE)
    hapsum.data = np.ones(hapsum.nnz)
    if alnmat.count is not None:
        cntmat = hapsum.transpose() * diags(alnmat.count, 0) * hapsum
    else:
        cntmat = hapsum.transpose() * hapsum
    return cntmat


def emase2db(sample_file, gene_file, db_file):
    start_time = time.time()

    LOG.info('Using sample file: {}'.format(sample_file))
    LOG.info('Using gene information file: {}'.format(gene_file))
    LOG.info('Generating db file: {}'.format(db_file))

    try:
        db_con = db_utils.init_db(db_file)
        chr_regex = re.compile('''(CHR(OMOSOME)?)?([1-9][0-9]*|X|Y|MT?)''', re.IGNORECASE)
        cursor = db_con.cursor()


        with tempfile.NamedTemporaryFile() as group_file_fd:
            group_file = group_file_fd.name

            temp_time = time.time()
            LOG.info("Inserting genes...")

            with open(gene_file, 'rb') as gene_fd:
                gene_info_tab_reader = csv.reader(gene_fd, delimiter='\t')
                for row in gene_info_tab_reader:
                    row = [x.strip() for x in row][:8]
                    gene_id, chromosome, start_pos_bp, end_pos_bp, strand, symbol, name, transcripts = row
                    chr_match = chr_regex.match(chromosome)
                    if chr_match:
                        chromosome = chr_match.group(3)
                        db_utils.add_gene_info(cursor, gene_id, chromosome, start_pos_bp, end_pos_bp, strand, symbol, name)
                        group_file_fd.write("{}\t{}\n".format(gene_id, "\t".join(transcripts.split(","))))

            LOG.info("Genes inserted in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                   utils.format_time(start_time, time.time())))

            temp_time = time.time()
            LOG.info("Inserting counts and edges...")

            with open(sample_file) as sample_fd:
                sample_reader = csv.reader(sample_fd, delimiter='\t')
                for row in sample_reader:
                    row = [x.strip() for x in row][:2]
                    sample_id, emase_file_name = row
                    db_utils.add_sample_id(cursor, sample_id, emase_file_name)

                    apm_temp_time = time.time()

                    LOG.info("Parsing: {}".format(emase_file_name))

                    alnmat = APM(h5file=emase_file_name, grpfile=group_file)
                    alnmat._bundle_inline(reset=True)
                    cntmat = get_num_shared_multireads(alnmat)

                    for nonzero_index_pair in np.transpose(cntmat.nonzero()):
                        i1, i2 = tuple(nonzero_index_pair)
                        read_count = cntmat[i1, i2]
                        if i1 == i2:
                            db_utils.add_gene_count_total(cursor, sample_id, alnmat.lname[i1], read_count)
                        else:
                            db_utils.add_gene_edge(cursor, sample_id, alnmat.lname[i1], alnmat.lname[i2], read_count)

                    LOG.info("Parsed in: {}, total time: {}".format(utils.format_time(apm_temp_time, time.time()),
                                                                    utils.format_time(start_time, time.time())))

            LOG.info("All counts and edges inserted in: {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                                   utils.format_time(start_time, time.time())))
        db_con.commit()

    except Exception as e:
        LOG.error('Error: {}'.format(str(e)))

    LOG.info("Database created in {}".format(utils.format_time(start_time, time.time())))


def parsefastqtest(input_file, chunks, dir):
    start_time = time.time()
    LOG.info('Using fastq file: {}'.format(input_file))
    barcode_utils.parse(input_file, chunks, dir)
    LOG.info("FAstq file parsed in {}".format(utils.format_time(start_time, time.time())))


def dumpec(ec_file):
    ECF = bin_file.ECFile(ec_file)


def ec2apm(ec_file, apm_file):
    bin_utils.bin2apm(ec_file, apm_file)
