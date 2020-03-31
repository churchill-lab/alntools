# -*- coding: utf-8 -*-
from collections import OrderedDict, namedtuple
from struct import pack
import multiprocessing
import os
import struct
import sys
import time

from six import iteritems

from Bio import bgzf
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix

import numpy as np
import pysam

from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from . import utils
from . import bin_utils

try:
    xrange
except NameError:
    xrange = range


LOG = utils.get_logger()


def parse_salmon_ec(salmon_dir, target_filename=None):

    LOG.debug('-------------------------------------------')
    LOG.debug('Parameters:')
    LOG.debug('  SALMON directory: {}'.format(salmon_dir))
    LOG.debug('  Target file: {}'.format(target_filename))
    LOG.debug('-------------------------------------------')

    if target_filename is not None:
        LOG.info('Reading in the reference transcript IDs from {}'.format(target_filename))
        transcripts_full = np.loadtxt(target_filename, dtype=str, delimiter='\t', usecols=(0,))

    transcript_idx = OrderedDict()
    haplotype_idx = OrderedDict()
    tidx2coord = OrderedDict()
    targets = list()
    quant_file = os.path.join(salmon_dir, 'quant.sf')
    salmon_file = os.path.join(salmon_dir, 'aux_info', 'eq_classes.txt')
    with open(salmon_file) as fh:
        LOG.info("Parsing {}".format(salmon_file))
        num_targets = int(fh.readline())
        num_ec = int(fh.readline())
        tidx = 0
        hidx = 0
        for _ in range(num_targets):
            curline = fh.readline()
            tkey = curline.rstrip()
            targets.append(tkey)
            tgt, hap = tkey.split('_')
            if not tgt in transcript_idx:
                transcript_idx[tgt] = tidx
                tidx += 1
            if not hap in haplotype_idx:
                haplotype_idx[hap] = hidx
                hidx += 1
        if target_filename is not None:
            for tgt in transcripts_full:
                if not tgt in transcript_idx:
                    transcript_idx[tgt] = tidx
                    tidx += 1
        for tidx, tkey in enumerate(targets):
            tgt, hap = tkey.split('_')
            tidx2coord[tidx] = [transcript_idx[tgt], haplotype_idx[hap]]
        num_haps = len(haplotype_idx)
        num_transcripts = len(transcript_idx)

        LOG.info('Reading in the effective transcript lengths from {}'.format(quant_file))
        targets = dict(zip(targets, np.arange(num_targets)))
        with open(quant_file) as qfh:
            qfh.readline()
            for curline in qfh:
                item = curline.rstrip().split('\t')
                tidx2coord[targets[item[0]]].append(float(item[2]))

        LOG.info('Creating EC alignment incidence matrix')
        rowid = list()
        colid = list()
        for h in range(num_haps):
            rowid.append(list())
            colid.append(list())

        ecidx = 0
        ec_counts = list()
        for curline in fh:
            item = curline.rstrip().split('\t')
            ec_counts.append(int(item[-1]))
            for ecitem in item[1:-1]:
                tidx, hidx, _ = tidx2coord[int(ecitem)]
                rowid[hidx].append(ecidx)
                colid[hidx].append(tidx)
            ecidx += 1

        data = list()
        for h in range(num_haps):
            data.append(coo_matrix((np.ones(len(rowid[h])), (rowid[h], colid[h])), shape=(num_ec, num_transcripts), dtype=int))
        for h in range(num_haps):
            data[h] = data[h].tocsr()

        alnmat = data[0]
        for h in range(1, num_haps):
           alnmat = alnmat + ((2 ** h) * data[h])

        LOG.info('Creating EC count matrix')
        ec_counts = np.array(ec_counts)
        cntmat = csr_matrix(ec_counts).transpose()  # Note: This is eventually a csc_matrix

        LOG.info('Creating objects of meta data')
        transcript_lengths = np.zeros((num_transcripts, num_haps), dtype=int)
        for _, tinfo in tidx2coord.items():
            transcript_lengths[tinfo[0], tinfo[1]] = tinfo[2]
        #transcript_lengths[transcript_lengths < 1] = 1.0
        haplotypes = np.array(list(haplotype_idx.keys()))
        transcripts = np.array(list(transcript_idx.keys()))

    LOG.info('Parsing complete')
    return transcripts, haplotypes, alnmat, cntmat, transcript_lengths


def convert(salmon_dir, ec_filename, sample, target_filename=None):

    LOG.debug('-------------------------------------------')
    LOG.debug('Parameters:')
    LOG.debug('  SALMON directory: {}'.format(salmon_dir))
    LOG.debug('  EC file: {}'.format(ec_filename))
    LOG.debug('  Sample: {}'.format(sample))
    LOG.debug('  Target file: {}'.format(target_filename))
    LOG.debug('-------------------------------------------')

    time0 = time.time()
    transcripts, haplotypes, alnmat, cntmat, transcript_lengths = parse_salmon_ec(salmon_dir, target_filename)
    LOG.info("{} parsed in {}".format(ec_filename, utils.format_time(time0, time.time())))

    time1 = time.time()
    LOG.info("Converting and storing to {}".format(ec_filename))
    bin_utils.save(ec_filename, sample, haplotypes, transcripts, transcript_lengths, alnmat, cntmat)
    LOG.info("{} created in {}".format(ec_filename, utils.format_time(time1, time.time())))
