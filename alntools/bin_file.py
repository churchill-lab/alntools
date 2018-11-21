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

    def get_ec_crs_dict(self):
        start_time = time.time()
        ecs = OrderedDict()
        for idx in xrange(len(self.n_matrix.indptr) - 1):
            ec_key = ','.join(map(str, self.n_matrix.getrow(idx).nonzero()[1]))
            #ec_key = self.a_matrix.getrow(idx).nonzero()[1]
            ecs[ec_key] = len(ecs)
        LOG.info("dict: {}".format(
                utils.format_time(start_time, time.time())))


    def toAPM(self):
        try:
            start_time = time.time()
            temp_time = time.time()

            num_haplotypes = len(self.haplotypes_idx)

            new_shape = (len(self.targets_idx),
                         num_haplotypes,
                         self.a_matrix.shape[0])

            LOG.debug('Shape={}'.format(new_shape))

            # final.ec.values -> the number of times this equivalence class has appeared

            ec_ids = [x for x in xrange(0, self.a_matrix.shape[0])]
            ec_arr = [[] for _ in xrange(0, num_haplotypes)]
            target_arr = [[] for _ in xrange(0, num_haplotypes)]

            for idx in xrange(self.a_matrix.shape[0]):
                # get the row values
                a_row = self.a_matrix.getrow(idx)

                # only need the columns (targets) that have values
                a_col = list(a_row.nonzero()[1])

                #print a_row
                #print a_col

                if num_haplotypes == 1:
                    # no need to decode
                    # print 'no decoding'
                    ec_arr[0].extend([idx] * a_row.nnz)
                    target_arr[0].extend(a_col)
                else:
                    #print 'decoding'
                    for x in a_col:
                        #print x
                        hap_values = utils.int_to_list(a_row[0, x], num_haplotypes)
                        #print hap_values

                        for i, h in enumerate(hap_values):
                            if h != 0:
                                ec_arr[i].append(idx)
                                target_arr[i].append(x)



            #for h in xrange(0, num_haplotypes):
            #    print h, len(ec_arr[h]), len(ec_ids), len(ECF.targets_idx)


            apm = APM(shape=new_shape,
                      haplotype_names=self.haplotypes_idx,
                      locus_names=self.targets_idx,
                      read_names=ec_ids,
                      sample_names=self.samples_idx)

            for h in xrange(0, num_haplotypes):
                d = np.ones(len(ec_arr[h]), dtype=np.int32)
                apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])), shape=(len(ec_ids), len(self.targets_idx)))

            apm.count = self.n_matrix

            LOG.info("APM Created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                utils.format_time(start_time, time.time())))

            return apm

        except KeyboardInterrupt as e:
            LOG.error("toAPM Error: {}".format(str(e)))
