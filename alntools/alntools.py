# -*- coding: utf-8 -*-

from collections import OrderedDict
import time
import os

from . import bam_utils
from . import utils

import numpy as np
from struct import pack

import emase
from emase import AlignmentPropertyMatrix as APM

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

    apm = emase.AlignmentPropertyMatrix(h5file=emase_filename)

    LOG.debug("EMASE APM created")

    em.target_list = list(apm.lname)
    em.target_dict = {target: idx for idx, target in enumerate(em.target_list)}

    em.haplotypes_list = apm.hname
    em.haplotypes_dict = {hap: idx for idx, hap in enumerate(em.haplotypes_list)}

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


class EC:
    def __init__(self, filename=None):
        self.version = -1
        self.filename = filename

        self.targets_list = []
        self.targets_dict = OrderedDict()

        self.haplotypes_list = []
        self.haplotypes_dict = OrderedDict()

        # Version 0

        self.reads_list = []
        self.reads_dict = OrderedDict()

        # Version 1

        self.ec_list = []
        self.ec_counts_list = []

        self.alignments = []


def parse_ec(file_in):
    """

    :param file_in:
    :return:
    """

    if not file_in:
        raise ValueError("empty file name, cannot load")

    f = open(file_in, 'rb')

    ec = EC(file_in)

    ec.version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

    if ec.version == 0:
        LOG.info("Version: 0, Reads")
    elif ec.version == 1:
        LOG.info("Version: 1, Equivalence Class")
    else:
        LOG.info("Unknown version, exiting")
        LOG.info("Exiting")
        return

    # TARGETS

    num_targets = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Target Count: {0:,}".format(num_targets))

    ec._targets_list = []
    ec._targets_dict = OrderedDict()

    for i in xrange(0, num_targets):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        target = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec.targets_dict[target] = i
        ec.targets_list.append(target)

    LOG.debug('ec.targets_list[0:10]={}'.format(str(ec.targets_list[0:10])))
    LOG.debug('len(ec.targets_list)={}'.format(len(ec.targets_list)))

    # HAPLOTYPES

    num_haplotypes = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Haplotype Count: {0:,}".format(num_haplotypes))

    ec.haplotypes_list = []
    ec.haplotypes_dict = OrderedDict()

    for i in xrange(0, num_haplotypes):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        haplotype = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec.haplotypes_dict[haplotype] = i
        ec.haplotypes_list.append(haplotype)

    if ec.version == 0:
        # READS

        ec.reads_list = []
        ec.reads_dict = OrderedDict()

        num_reads = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Read Count: {0:,}".format(num_reads))

        for i in xrange(0, num_reads):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            read_id = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            ec.reads_dict[read_id] = i
            ec.reads_list.append(read_id)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        ec.alignments = []

        for i in xrange(0, num_alignments*3, 3):
            read_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec.alignments.append((read_index, target_index, bit_flag))
    else:
        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Equivalance Class Count: {0:,}".format(num_ec))

        ec.ec_list = [x for x in xrange(0, num_ec)]
        ec.ec_counts_list = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        ec.alignments = []

        for i in xrange(0, num_alignments*3, 3):
            ec_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec.alignments.append((ec_index, target_index, bit_flag))

    return ec


def ec2emase(file_in, file_out):
    """

    :param file_in:
    :param file_out:
    :return:
    """
    start_time = time.time()

    LOG.info('Parsing EC file: {}'.format(file_in))
    temp_time = time.time()
    ec = parse_ec(file_in)

    LOG.info("EC parsed in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                        utils.format_time(start_time, time.time())))

    new_shape = (len(ec.targets_list), len(ec.haplotypes_list), len(ec.ec_list))

    LOG.info('Constructing APM structure...')
    LOG.debug('Shape={}'.format(new_shape))

    apm = APM(shape=new_shape, haplotype_names=ec.haplotypes_list, locus_names=ec.targets_list, read_names=ec.ec_list)

    LOG.debug('ec.haplotypes_list={}'.format(str(ec.haplotypes_list)))
    LOG.debug('ec.targets_list[0:10]={}'.format(str(ec.targets_list[0:10])))
    LOG.debug('ec.ec_list[0:10]={}'.format(str(ec.ec_list[0:10])))

    temp_time = time.time()

    # counts -> the number of times this equivalence class has appeared
    apm.count = ec.ec_counts_list

    num_haplotypes = len(ec.haplotypes_list)

    try:

        for alignment in ec.alignments:
            # LOG.debug(str(alignment))
            ec_index = alignment[0]
            target_index = alignment[1]
            temp_bits = alignment[2]

            if temp_bits == 0:
                continue

            bits = utils.int_to_list(temp_bits, num_haplotypes)
            for i, bit in enumerate(bits):
                if bit:
                    # lid, hid, rid, value
                    apm.set_value(target_index, i, ec_index, 1)
    except Exception as e:
        LOG.error('Error: {}'.format(str(e)))

    LOG.info("APM Created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                        utils.format_time(start_time, time.time())))

    temp_time = time.time()
    LOG.info("Flushing to disk...")
    apm.finalize()
    apm.save(file_out, title='bam2ec')
    LOG.info("{} created in {}, total time: {}".format(file_out,
                                                   utils.format_time(temp_time, time.time()),
                                                   utils.format_time(start_time, time.time())))


def bam2ec(bam_filename, ec_filename, num_chunks=0, target_filename=None, temp_dir=None, range_filename=None):
    bam_utils.convert(bam_filename, ec_filename, num_chunks=num_chunks, target_filename=target_filename, emase=False, temp_dir=temp_dir, range_filename=range_filename)


def bam2emase(bam_filename, emase_filename, num_chunks=0, target_filename=None, temp_dir=None):
    bam_utils.convert(bam_filename, emase_filename, num_chunks=num_chunks, target_filename=target_filename, emase=True, temp_dir=temp_dir)


def split_bam(bam_filename, num_chunks, directory=None):
    bam_utils.split_bam(bam_filename, num_chunks, directory)


def generate_bam_ranges(input_files, range_filename, target_filename=None, temp_dir=None):
    bam_utils.generate_bam_ranges(input_files, range_filename, target_filename, temp_dir)
