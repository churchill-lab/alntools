# -*- coding: utf-8 -*-
from collections import OrderedDict
from struct import pack, unpack, calcsize
from ntpath import basename

import csv
import gzip
import os
import re
import tempfile
import time

import numpy as np

from scipy.sparse import coo_matrix, diags, csr_matrix, csc_matrix

from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from . import bam_utils
from . import db_utils
from . import bam_utils_multisample
from . import barcode_utils
from . import utils

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


class EMASE:
    def __init__(self, filename=None):
        self.filename = filename



        self.target_list = []
        self.target_dict = OrderedDict()

        self.haplotypes_list = []
        self.haplotypes_dict = OrderedDict()

        # Version 0

        self.reads_list = []
        self.reads_dict = OrderedDict()

        # Version 1

        self.ec_list = []
        self.ec_counts_list = []

        self.alignments = []


def parse_emase(emase_filename):
    """
    Parse EMASE file.

    :param emase_filename: Name of the EMASE file
    :return: EMASE class
    """
    if not emase_filename:
        raise ValueError("empty file name, cannot load")

    em = EMASE(emase_filename)

    LOG.debug("Creating EMASE APM")

    apm = APM(h5file=emase_filename)

    LOG.debug("EMASE APM created")

    em.target_list = list(apm.lname)
    em.target_dict = {target: idx for idx, target in enumerate(em.target_list)}

    em.haplotypes_list = apm.hname
    em.haplotypes_dict = {hap: idx for idx, hap in enumerate(em.haplotypes_list)}

    if apm.rname is None:
        em.ec_list = np.arange(apm.num_reads)
    else:
        em.ec_list = list(apm.rname)
    em.ec_counts_list = list(apm.count)

    LOG.debug('Adding {:,} elements'.format(len(em.ec_list)))

    for ec_idx, ec in enumerate(em.ec_list):
        for target, target_idx in em.target_dict.iteritems():
            #LOG.debug('ec_idx = {}, target = {}'.format(ec_idx, target))
            bits = []

            for hap, hap_idx in em.haplotypes_dict.iteritems():
                bits.append(apm.data[hap_idx][ec_idx, target_idx])


            em.alignments.append((ec_idx, target_idx, utils.list_to_int(bits)))

    LOG.debug("EMASE file parsed")

    return em


def emase2ec(emase_filename, ec_filename):
    """
    Convert the EMASE file to an EC file.

    :param emase_filename: name of the EMASE file
    :param ec_filename: name of the EC file, version 1
    :return:
    """
    start_time = time.time()

    LOG.info('Parsing Emase file: {}'.format(emase_filename))
    temp_time = time.time()
    emase_data = parse_emase(emase_filename)

    LOG.info("Emase file parsed in {}, total time: {}".format(emase_filename,
                                                              utils.format_time(temp_time, time.time()),
                                                              utils.format_time(start_time, time.time())))

    try:
        LOG.info("Generating EC file...")

        f = open(ec_filename, "wb")

        # version
        f.write(pack('<i', 1))
        LOG.debug("1\t# VERSION")

        # targets
        LOG.info("{:,}\t# NUMBER OF TARGETS".format(len(emase_data.target_list)))
        f.write(pack('<i', len(emase_data.target_list)))
        for main_target, idx in emase_data.target_dict.iteritems():
            LOG.debug("{:,}\t{}\t# {:,}".format(len(main_target), main_target, idx))
            f.write(pack('<i', len(main_target)))
            f.write(pack('<{}s'.format(len(main_target)), main_target))

        # haplotypes
        LOG.info("{:,}\t# NUMBER OF HAPLOTYPES".format(len(emase_data.haplotypes_list)))
        f.write(pack('<i', len(emase_data.haplotypes_list)))
        for idx, hap in enumerate(emase_data.haplotypes_list):
            LOG.debug("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
            f.write(pack('<i', len(hap)))
            f.write(pack('<{}s'.format(len(hap)), hap))

        # equivalence classes
        LOG.info("{:,}\t# NUMBER OF EQUIVALANCE CLASSES".format(len(emase_data.ec_list)))
        f.write(pack('<i', len(emase_data.ec_list)))
        for idx, k in enumerate(emase_data.ec_counts_list):
            # k is the count
            LOG.debug("{:,}\t# {:,}".format(k, idx))
            f.write(pack('<i', k))

        LOG.info("Determining mappings...")

        # equivalence class mappings
        LOG.info("{:,}\t# NUMBER OF EQUIVALANCE CLASS MAPPINGS".format(len(emase_data.alignments)))
        f.write(pack('<i', len(emase_data.alignments)))

        for idx, alignment in enumerate(emase_data.alignments):
            LOG.debug("{}\t{}\t{}\t# {}\t{}".format(alignment[0], alignment[1], alignment[2], emase_data.target_list[alignment[1]], utils.int_to_list(alignment[2], len(emase_data.haplotypes_list))))
            f.write(pack('<i', alignment[0]))
            f.write(pack('<i', alignment[1]))
            f.write(pack('<i', alignment[2]))

        f.close()
    except Exception as e:
        LOG.error("Error: {}".format(str(e)))


class ECFile:
    def __init__(self, filename=None):
        self.filename = filename

        # TODO: implement
        self.header = None
        self.headers = None

        self.format = -1

        # key   = haplotype name
        # value = haplotype insertion order
        self.haplotypes = OrderedDict()
        self.haplotypes_idx = []

        # key   = target name
        # value = target insertion order
        self.targets = OrderedDict()
        self.targets_idx = []
        self.targets_lengths = OrderedDict()

        # key   = sample (cr) name
        # value = sample (cr) insertion order
        self.samples = OrderedDict()
        self.samples_idx = []


        #
        # Format 0 (NO LONGER SUPPORTED)

        #
        # Format 1 will be stored in Format 2
        #

        #
        # SECTION: "A" Matrix
        #
        # "A" Matrix format is EC (rows) by Transcripts (columns) with
        # each value being the HAPLOTYPE flag.
        #
        # Instead of storing a "dense" matrix, we store a "sparse"
        # matrix utilizing Compressed Sparse Row (CSR) format.
        #
        # NOTE:
        #     HAPLOTYPE flag is an integer that denotes which haplotype
        #     (allele) a read aligns to given an EC. For example, 00,
        #     01, 10, and 11 can specify whether a read aligns to the
        #     1st and/or 2nd haplotype of a transcript.  These binary
        #     numbers are converted to integers - 0, 1, 2, 3 - and
        #     stored as the haplotype flag.
        #

        self.a_matrix = None

        #
        # SECTION: "N" Matrix
        #
        # "N" Matrix format is EC (rows) by CRS (columns) with
        # each value being the EC count.
        #
        # Instead of storing a "dense" matrix, we store a "sparse"
        # matrix utilizing Compressed Sparse Column (CSC) format.
        #

        self.n_matrix = None

        self.__load__()

    def __load__(self):
        if not self.filename:
            raise ValueError("empty file name, cannot load")

        LOG.debug("Parsing: {}".format(self.filename))

        _i = calcsize("<i")
        _s = calcsize("<s")

        start_time = time.time()

        # attempt to see if this is gzipped (version 2)
        with gzip.open(self.filename, 'rb') as f:

            self.format = unpack('<i', f.read(_i))[0]

            if self.format == 0:
                LOG.error("Version: 0, Reads")
                LOG.error('Version no longer supported')
                raise ValueError("Unsupported Version 0")
            elif self.format == 1:
                LOG.error("Version: 1, Equivalence Class")
            elif self.format == 2:
                LOG.error("Version: 2, Multisample")
            else:
                LOG.error("Unknown version: {}, exiting".format(self.format))
                LOG.error("Exiting")
                return

            #
            # SECTION: HAPLOTYPES
            #     [# of HAPLOTYPES = H]
            #     [H1 text length][H1 text]
            #     ...
            #     [HH text length][HH text]
            #
            # Example:
            #     8
            #     1 A
            #     1 B
            #     1 C
            #     1 D
            #     1 E
            #     1 F
            #     1 G
            #     1 H
            #

            temp_time = time.time()

            self.haplotypes = OrderedDict()
            self.haplotypes_idx = []

            num_haplotypes = unpack('<i', f.read(_i))[0]
            LOG.error("Haplotype Count: {0:,}".format(num_haplotypes))

            for i in xrange(0, num_haplotypes):
                str_len = unpack('<i', f.read(_i))[0]
                haplotype = unpack('<{}s'.format(str_len), f.read(_s * str_len))[0]
                self.haplotypes[haplotype] = i
                self.haplotypes_idx.append(haplotype)
                LOG.debug("{} {}".format(i, haplotype))

            LOG.info("Haplotypes extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

            #
            # SECTION: TARGETS
            #     [# of TARGETS = T]
            #     [T1 text length][T1 text][HAP 1 length] ... [HAP H length]
            #     ...
            #     [TT text length][TT text][HAP 1 length] ... [HAP H length]
            #
            # Example:
            #     80000
            #     18 ENSMUST00000156068 234 235
            #     18 ENSMUST00000209341 1054 1054
            #     ...
            #     18 ENSMUST00000778019 1900 1899
            #

            temp_time = time.time()

            self.targets = OrderedDict()
            self.targets_idx = []
            self.targets_lengths = OrderedDict()

            num_targets = unpack('<i', f.read(_i))[0]
            LOG.error("Target Count: {0:,}".format(num_targets))

            for i in xrange(0, num_targets):
                str_len = unpack('<i', f.read(_i))[0]
                target = unpack('<{}s'.format(str_len), f.read(_s * str_len))[0]
                self.targets[target] = i
                self.targets_idx.append(target)
                LOG.debug("{} {}".format(i, target))

                hap_length = OrderedDict()
                for haplotype, idx in self.haplotypes.iteritems():
                    hap_length[haplotype] = unpack('<i', f.read(_i))[0]

                self.targets_lengths[target] = hap_length

            LOG.info("Targets extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

            if self.format == 1:
                #
                # SECTION: "A" Matrix
                #
                # "A" Matrix format is EC (rows) by Transcripts (columns) with
                # each value being the HAPLOTYPE flag.
                #
                # Instead of storing a "dense" matrix, we store a "sparse"
                # matrix utilizing Compressed Sparse Row (CSR) format.
                #
                # NOTE:
                #     HAPLOTYPE flag is an integer that denotes which haplotype
                #     (allele) a read aligns to given an EC. For example, 00,
                #     01, 10, and 11 can specify whether a read aligns to the
                #     1st and/or 2nd haplotype of a transcript.  These binary
                #     numbers are converted to integers - 0, 1, 2, 3 - and
                #     stored as the haplotype flag.
                #

                indptr_length = unpack('<i', f.read(_i))[0]
                nnz = unpack('<i', f.read(_i))[0]

                LOG.debug("A MATRIX INDPTR Length: {0:,}".format(indptr_length))
                LOG.debug("A MATRIX NNZ: {0:,}".format(nnz))

                indptr = np.array(unpack('<{}i'.format(indptr_length), f.read(_i * indptr_length)), dtype=np.int32)
                indices = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)
                data = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)

                self.a_matrix = csr_matrix((data, indices, indptr))

                #
                # EC
                #

                num_ec = unpack('<i', f.read(_i))[0]
                print("EC Count: {0:,}".format(num_ec))

                ec = list(unpack('<{}i'.format(num_ec), f.read(_i * num_ec)))

                # convert to Format 2, "N" Matrix
                self.samples = OrderedDict()
                self.samples[basename(self.filename)] = len(self.samples)
                self.samples_idx.append(basename(self.filename))

                self.n_matrix = csc_matrix(ec)

            elif self.format == 2:
                #
                # CR
                #

                temp_time = time.time()

                self.samples = OrderedDict()

                num_samples = unpack('<i', f.read(_i))[0]
                LOG.error("Sample Count: {0:,}".format(num_samples))

                for i in xrange(0, num_samples):
                    str_len = unpack('<i', f.read(_i))[0]
                    sample = unpack('<{}s'.format(str_len), f.read(_s * str_len))[0]
                    self.samples[sample] = i
                    self.samples_idx.append(sample)
                    LOG.debug("{} {}".format(i, sample))

                LOG.info("Samples extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                #
                # SECTION: "A" Matrix
                #
                # "A" Matrix format is EC (rows) by Transcripts (columns) with
                # each value being the HAPLOTYPE flag.
                #
                # Instead of storing a "dense" matrix, we store a "sparse"
                # matrix utilizing Compressed Sparse Row (CSR) format.
                #
                # NOTE:
                #     HAPLOTYPE flag is an integer that denotes which haplotype
                #     (allele) a read aligns to given an EC. For example, 00,
                #     01, 10, and 11 can specify whether a read aligns to the
                #     1st and/or 2nd haplotype of a transcript.  These binary
                #     numbers are converted to integers - 0, 1, 2, 3 - and
                #     stored as the haplotype flag.
                #

                temp_time = time.time()
                a_time = time.time()

                indptr_length = unpack('<i', f.read(_i))[0]
                nnz = unpack('<i', f.read(_i))[0]

                LOG.error("A MATRIX INDPTR Length: {0:,}".format(indptr_length))
                LOG.error("A MATRIX NNZ: {0:,}".format(nnz))

                indptr = np.array(unpack('<{}i'.format(indptr_length), f.read(_i * indptr_length)), dtype=np.int32)

                LOG.info("indptr extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                indices = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)

                LOG.info("indices extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                data = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)

                LOG.info("data extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                self.a_matrix = csr_matrix((data, indices, indptr))

                LOG.info("A matrix created in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                LOG.info("A matrix extraction and creation: {}".format(
                        utils.format_time(a_time, time.time())))

                #
                # SECTION: "N" Matrix
                #
                # "N" Matrix format is EC (rows) by CRS (columns) with
                # each value being the EC count.
                #
                # Instead of storing a "dense" matrix, we store a "sparse"
                # matrix utilizing Compressed Sparse Column (CSC) format.
                #
                # NOTE: Since this is CSC instead of CSR
                #
                # CSR
                # - indptr points to row starts in indices and data
                # - indices is array of column indices
                # - data is array of corresponding nonzero values
                #
                # CSC
                # - indptr points to column starts in indices and data
                # - indices is array of row indices
                # - data is array of corresponding nonzero values
                #

                temp_time = time.time()
                n_time = time.time()

                indptr_length = unpack('<i', f.read(_i))[0]
                nnz = unpack('<i', f.read(_i))[0]

                LOG.error("N MATRIX INDPTR Length: {0:,}".format(indptr_length))
                LOG.error("N MATRIX NNZ: {0:,}".format(nnz))

                indptr = np.array(unpack('<{}i'.format(indptr_length), f.read(_i * indptr_length)), dtype=np.int32)

                LOG.info("indptr extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                indices = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)

                LOG.info("indices extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                data = np.array(unpack('<{}i'.format(nnz), f.read(_i * nnz)), dtype=np.int32)

                LOG.info("data extracted in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                temp_time = time.time()

                self.n_matrix = csc_matrix((data, indices, indptr))

                LOG.info("N matrix created in {}, total time: {}".format(
                        utils.format_time(temp_time, time.time()),
                        utils.format_time(start_time, time.time())))

                LOG.info("N matrix extraction and creation: {}".format(
                        utils.format_time(n_time, time.time())))

                LOG.info("All data parsed in: {}".format(
                        utils.format_time(start_time, time.time())))

    def get_samples(self):
        pass

    def get_ec_dict(self):
        start_time = time.time()
        ecs = OrderedDict()
        for idx in xrange(len(self.a_matrix.indptr) - 1):
            ec_key = ','.join(map(str, self.a_matrix.getrow(idx).nonzero()[1]))
            #ec_key = self.a_matrix.getrow(idx).nonzero()[1]
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))

    def get_ec_dict2(self):
        start_time = time.time()
        ecs = {}
        for idx in xrange(len(self.a_matrix.indptr) - 1):
            ec_key = ','.join(map(str, self.a_matrix.getrow(idx).nonzero()[1]))
            #ec_key = self.a_matrix.getrow(idx).nonzero()[1]
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))

    def get_ec_dict3(self):
        start_time = time.time()
        ecs = OrderedDict()
        for idx in xrange(len(self.a_matrix.indptr) - 1):
            #ec_key = ','.join(map(str, self.a_matrix.getrow(idx).nonzero()[1]))
            ec_key = tuple(self.a_matrix.getrow(idx).nonzero()[1].tolist())
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))

    def get_ec_dict4(self):
        start_time = time.time()
        ecs = {}
        for idx in xrange(len(self.a_matrix.indptr) - 1):
            #ec_key = ','.join(map(str, self.a_matrix.getrow(idx).nonzero()[1]))
            ec_key = tuple(self.a_matrix.getrow(idx).nonzero()[1].tolist())
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))



    def get_ec_crs_dict(self):
        start_time = time.time()
        ecs = OrderedDict()
        for idx in xrange(len(self.n_matrix.indptr) - 1):
            ec_key = ','.join(map(str, self.n_matrix.getrow(idx).nonzero()[1]))
            #ec_key = self.a_matrix.getrow(idx).nonzero()[1]
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))
