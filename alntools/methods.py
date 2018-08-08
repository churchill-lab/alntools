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


class EC:
    def __init__(self, filename=None):
        self.version = -1
        self.filename = filename

        self.haplotypes = OrderedDict()
        self.targets = OrderedDict()

        #
        # Version 0
        #

        self.reads_list = []
        self.reads_dict = OrderedDict()

        #
        # Version 1
        #

        # "A" Matrix (Compressed Sparse Row format)

        self.a_indptr = None
        self.a_indices = None
        self.a_data = None

        # Equivalence Classes
        self.ec = []

        #
        # Version 2
        #

        # "A" Matrix (Compressed Sparse Row format)

        self.ec_list = []
        self.ec_counts_list = []
        self.alignments = []

def parse_ec(file_in, detail=False):
    """

    :param file_in:
    :return:
    """

    if not file_in:
        raise ValueError("empty file name, cannot load")

    ec = EC(file_in)

    # attempt to see if this is gzipped (version 2)
    try:
        f = gzip.open(file_in)
        ec.version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    except Exception as e:
        print('error={}'.format(str(e)))
        f = open(file_in, 'rb')
        ec.version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

    LOG.debug("Parsing: {}".format(file_in))

    if ec.version == 0:
        LOG.error("Version: 0, Reads")
    elif ec.version == 1:
        LOG.error("Version: 1, Equivalence Class")
    elif ec.version == 2:
        LOG.error("Version: 2, Multisample")
    else:
        LOG.error("Unknown version: {}, exiting".format(ec.version))
        LOG.error("Exiting")
        return

    # TARGETS

    num_targets = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.error("Target Count: {0:,}".format(num_targets))

    for i in xrange(0, num_targets):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        target = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec.targets_dict[target] = i
        ec.targets_list.append(target)
        LOG.info("{} {}".format(i, target))

    #LOG.debug('ec.targets_list[0:10]={}'.format(str(ec.targets_list[0:10])))
    #LOG.debug('len(ec.targets_list)={}'.format(len(ec.targets_list)))

    # HAPLOTYPES

    num_haplotypes = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.error("Haplotype Count: {0:,}".format(num_haplotypes))

    ec.haplotypes_list = []
    ec.haplotypes_dict = OrderedDict()

    for i in xrange(0, num_haplotypes):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        haplotype = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec.haplotypes_dict[haplotype] = i
        ec.haplotypes_list.append(haplotype)
        LOG.info("{} {}".format(i, haplotype))

    if ec.version == 0:
        # READS

        ec.reads_list = []
        ec.reads_dict = OrderedDict()

        num_reads = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Read Count: {0:,}".format(num_reads))

        for i in xrange(0, num_reads):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            read_id = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            ec.reads_dict[read_id] = i
            ec.reads_list.append(read_id)
            LOG.info("{} {}".format(i, read_id))

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        ec.alignments = []

        for i in xrange(0, num_alignments*3, 3):
            read_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec.alignments.append((read_index, target_index, bit_flag))
            LOG.info("{} {} {}  # {} {} {} ".format(read_index, target_index, bit_flag, ec.reads_list[read_index], ec.targets_list[target_index], bit_flag))
    elif ec.version == 1:

        indptr_length = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        nnz = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]




        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Equivalence Class Count: {0:,}".format(num_ec))

        ec.ec_list = [x for x in xrange(0, num_ec)]
        ec.ec_counts_list = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        ec.alignments = []

        for i in xrange(0, num_alignments*3, 3):
            ec_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec.alignments.append((ec_index, target_index, bit_flag))
            LOG.info("{} {} {}  # {} {} {} ".format(ec_index, target_index, bit_flag, ec.ec_list[ec_index], ec.targets_list[target_index], bit_flag))

    elif ec.version == 2:
        #
        # CR
        #

        cr_dict = OrderedDict()
        cr_list = []

        num_crs = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("CR Count: {0:,}".format(num_crs))

        for i in xrange(0, num_crs):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            cr = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            cr_dict[cr] = i
            cr_list.append(cr)
            LOG.info("{} {}".format(i, cr))

        # "N" MATRIX
        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Equivalence Class Count: {0:,}".format(num_ec))

        num_nnz = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Non-zero Count: {0:,}".format(num_nnz))

        # row offsets
        indptr = np.fromfile(f, dtype=np.dtype('i'), count=num_ec + 1)
        LOG.error(indptr)

        # columns
        indices = np.fromfile(f, dtype=np.dtype('i'), count=num_nnz)
        LOG.error(indices)

        # data
        data = np.fromfile(f, dtype=np.dtype('i'), count=num_nnz)
        LOG.error(data)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.error("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        LOG.error(temp_alignments)
        ec.alignments = []

        for i in xrange(0, num_alignments*3, 3):
            ec_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec.alignments.append((ec_index, target_index, bit_flag))
            LOG.info("{} {} {}".format(ec_index, target_index, bit_flag))

        LOG.error('done version 3')

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

    #LOG.debug('ec.haplotypes_list={}'.format(str(ec.haplotypes_list)))
    #LOG.debug('ec.targets_list[0:10]={}'.format(str(ec.targets_list[0:10])))
    #LOG.debug('ec.ec_list[0:10]={}'.format(str(ec.ec_list[0:10])))

    temp_time = time.time()

    # counts -> the number of times this equivalence class has appeared
    apm.count = ec.ec_counts_list

    num_haplotypes = len(ec.haplotypes_list)
    ec_arr = [[] for _ in xrange(0, num_haplotypes)]
    target_arr = [[] for _ in xrange(0, num_haplotypes)]

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
                    #apm.set_value(target_index, i, ec_index, 1)
                    ec_arr[i].append(ec_index)
                    target_arr[i].append(target_index)

        for h in xrange(0, num_haplotypes):
            d = np.ones(len(ec_arr[h]))
            apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])), shape=(len(ec.ec_list), len(ec.targets_list)))

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


def dumpec(ec_filename, detail=False):
    parse_ec(ec_filename, detail)


def bam2ec(bam_filename, ec_filename, chunks=0, directory=None, number_processes=-1, range_filename=None):
    bam_utils.convert(bam_filename, ec_filename, None, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2emase(bam_filename, emase_filename, chunks=0, directory=None, number_processes=-1, range_filename=None):
    bam_utils.convert(bam_filename, None, emase_filename, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2both(bam_filename, ec_filename, emase_filename, chunks=0, directory=None, number_processes=-1, range_filename=None):
    bam_utils.convert(bam_filename, ec_filename, emase_filename, num_chunks=chunks, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2ec_multisample(bam_filename, ec_filename, chunks=0, minimum_count=-1, directory=None, number_processes=-1, range_filename=None):
    bam_utils_multisample.convert(bam_filename, ec_filename, None, num_chunks=chunks, minimum_count=minimum_count, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


def bam2both_multisample(bam_filename, emase_filename, chunks=0, minimum_count=-1, directory=None, number_processes=-1, range_filename=None):
    bam_utils_multisample.convert(bam_filename, None, emase_filename, num_chunks=chunks, minimum_count=minimum_count, number_processes=number_processes, temp_dir=directory, range_filename=range_filename)


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
