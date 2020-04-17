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
from scipy.sparse import coo_matrix, csc_matrix

import numpy as np
import pysam

from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from . import utils

try:
    xrange
except NameError:
    xrange = range


LOG = utils.get_logger()
BAM_HEADER = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00'
BAM_EOF = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'

parse_fields = ["header_size", "begin_read_offset", "begin_read_size", "file_offset", "file_bytes", "end_read_offset",
                "end_read_size"]
ParseRecord = namedtuple("ParseRecord", parse_fields)


class ConvertParams(object):
    slots = ['input_file', 'temp_dir', 'process_id', 'track_ranges', 'data']

    def __init__(self):
        self.input_file = None
        self.temp_dir = None
        self.process_id = None
        self.track_ranges = False

        # tuple of (idx, ParseRecord)
        self.data = []

    def __str__(self):
        return "Input: {}\nProcess ID: {}\nData: {}".format(self.input_file, self.process_id, self.data)


class ConvertResults(object):
    slots = ['valid_alignments', 'all_alignments', 'ec', 'unique_reads', 'init', 'tid_ranges']

    def __init__(self):
        self.valid_alignments = None
        self.all_alignments = None
        self.ec = None
        self.unique_reads = None
        self.init = False
        self.tid_ranges = None


class RangeParams(object):
    slots = ['input_file', 'temp_dir', 'process_id']

    def __init__(self):
        self.input_file = None
        self.temp_dir = None
        self.process_id = None

    def __str__(self):
        return "Input: {}\nProcess ID: {}".format(self.input_file, self.process_id)


class RangeResults(object):
    slots = ['main_targets', 'haplotypes', 'init', 'tid_ranges']

    def __init__(self):
        self.main_targets = None
        self.haplotypes = None
        self.init = False
        self.tid_ranges = None


def get_header_size(bam_filename):
    """

    :param bam_filename:
    :return:
    """
    #
    # grab header
    #
    alignment_file = pysam.AlignmentFile(bam_filename)
    header_size = alignment_file.tell() >> 16
    alignment_file.close()
    return header_size


def fix_bam(filename):
    """
    Make sure the EOF marker is present.

    :param filename: the name of the BAME file
    :return: Nothing
    """
    if not os.path.isfile(filename):
        sys.exit("Missing file {}".format(filename))

    size = os.path.getsize(filename)
    h = open(filename, "rb")

    # Check it looks like a BGZF file
    # (could still be GZIP'd, in which case the extra block is harmless)
    data = h.read(len(BAM_HEADER))

    if data != BAM_HEADER:
        raise Exception("File {} is not a BAM file".format(filename))

    # Check if it has the EOF already
    h.seek(size - 28)
    data = h.read(28)
    h.close()

    if data != BAM_EOF:
        # Adding EOF block
        h = open(filename, "ab")
        h.write(BAM_EOF)
        h.close()


def validate_bam(filename):
    if not os.path.isfile(filename):
        sys.exit("Missing file {}".format(filename))

    size = os.path.getsize(filename)
    h = open(filename, "rb")

    # Check it looks like a BGZF file
    # (could still be GZIP'd, in which case the extra block is harmless)
    data = h.read(len(BAM_HEADER))

    if data != BAM_HEADER:
        raise Exception("File {} is not a BAM file".format(filename))

    # Check if it has the EOF already
    h.seek(size - 28)
    data = h.read(28)
    h.close()

    if data != BAM_EOF:
        raise Exception("File {} has bad EOF".format(filename))


def chunk_bam_file(bam_filename, new_filename, parse_rec):
    """
    Create a new BAM file from an existing one.

    :param str bam_filename: the name of the original BAM file
    :param str new_filename: the name of the new BAM file
    :param class:`ParseRecord` parse_rec: the information containing where to extract
    :return:
    """
    try:
        os.remove(new_filename)
    except Exception as e:
        pass

    # copy the header from original BAM file to new
    utils.bytes_from_file(bam_filename, new_filename, 0, parse_rec.header_size)

    if parse_rec.begin_read_offset > 0:
        # if there are reads before a chunk offset, we need to extract them
        b = bgzf.BgzfReader(bam_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.begin_read_offset)
        b2.write(b.read(parse_rec.begin_read_size))
        b2.close()
        truncate_bam_file(new_filename)

    # grab bgzf chunks from the OLD BAM file and append to NEW BAM file
    bytes_from_file_bam(bam_filename, new_filename, parse_rec.file_offset, parse_rec.file_bytes)

    if parse_rec.end_read_offset > 0:
        # if there are reads after a chunk offset, we need to extract them
        b = bgzf.BgzfReader(bam_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.end_read_offset)
        b2.write(b.read(parse_rec.end_read_size))
        b2.close()

    # fix the bam EOF if needed
    fix_bam(new_filename)


def process_convert_bam(cp):
    """

    :return:
    """
    LOG.debug('Process ID: {}, Input File: {}'.format(cp.process_id, cp.input_file))

    validate_bam(cp.input_file)

    LOG.debug('File Valid: {}, Input File: {}'.format(cp.process_id, cp.input_file))



    # reference_id: the reference sequence number as defined in the header
    # reference_name: name (None if no AlignmentFile is associated)

    # ec = equivalence class
    #      the KEY is a comma separated string of reference_ids
    #      the VALUE is the number of times this equivalence class has appeared
    ec = OrderedDict()

    # ec_idx = lookup to ec
    #          the KEY is a comma separated string of reference_ids
    #          the VALUE is a number specifying the insertion order of
    #               the KEY value in ec

    # unique reads
    unique_reads = {}

    # times encountering new read id
    read_id_switch_counter = 0
    same_read_target_counter = 0

    ranges = {}

    all_alignments = 0
    valid_alignments = 0
    temp_name = os.path.join(cp.temp_dir, '_bam2ec.')

    try:
        for file_info_data in cp.data:
            reference_id = None
            reference_ids = []

            try:
                idx = file_info_data[0]
                parse_record = file_info_data[1]

                # must create the file
                temp_file = "{}{}.bam".format(temp_name, idx)
                LOG.debug("Process ID: {}, Creating alignment file: {}".format(cp.process_id, temp_file))
                utils.delete_file(temp_file)
                chunk_bam_file(cp.input_file, temp_file, parse_record)
                LOG.debug("Process ID: {}, Opening alignment file: {}".format(cp.process_id, temp_file))

                alignment_file = pysam.AlignmentFile(temp_file)

                # query template name
                query_name = None

                while True:
                    alignment = alignment_file.__next__()

                    all_alignments += 1

                    # if alignment.flag == 4 or alignment.is_unmapped:
                    if alignment.is_unmapped:
                        continue

                    # added for paired end files
                    if alignment.is_paired:
                        if alignment.is_read2 or not alignment.is_proper_pair or alignment.reference_id != alignment.next_reference_id or alignment.next_reference_start < 0:
                            continue

                    valid_alignments += 1

                    # reference_sequence_name = Column 3 from file, the Reference NAME (EnsemblID_Haplotype)
                    # reference_id = the reference_id, which is 0 or a positive integer mapping to entries
                    #       within the sequence dictionary in the header section of a BAM file
                    # main_target = the Ensembl id of the transcript

                    reference_sequence_name = alignment.reference_name
                    reference_id = str(alignment.reference_id)

                    if cp.track_ranges:
                        min_max = ranges.get(reference_id, (100000000000, -1))
                        n = min(min_max[0], alignment.reference_start)
                        x = max(min_max[1], alignment.reference_start)
                        ranges[reference_id] = (n, x)

                    # query_name = Column 1 from file, the Query template NAME
                    if query_name is None:
                        query_name = alignment.query_name

                        i = query_name.find(' ')
                        if i > 0:
                            query_name = query_name[:i]

                    try:
                        unique_reads[query_name] += 1
                    except KeyError:
                        unique_reads[query_name] = 1

                    alignment_query_name = alignment.query_name
                    i = alignment_query_name.find(' ')
                    if i > 0:
                        alignment_query_name = alignment_query_name[:i]

                    if query_name != alignment_query_name:
                        ec_key = ','.join(sorted(reference_ids))

                        try:
                            ec[ec_key] += 1
                        except KeyError:
                            ec[ec_key] = 1

                        query_name = alignment.query_name
                        i = query_name.find(' ')
                        if i > 0:
                            query_name = query_name[:i]

                        reference_ids = [reference_id]
                        read_id_switch_counter += 1
                    else:
                        if reference_id not in reference_ids:
                            reference_ids.append(reference_id)
                        else:
                            same_read_target_counter += 1

                    if all_alignments % 100000 == 0:
                        LOG.debug("Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}, with {:,} equivalence classes".format(cp.process_id, temp_file, valid_alignments, all_alignments, len(ec)))

            except StopIteration:
                LOG.info(
                    "DONE Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}, with {:,} equivalence classes".format(
                        cp.process_id, temp_file, valid_alignments, all_alignments, len(ec)))

            #LOG.info("{0:,} alignments processed, with {1:,} equivalence classes".format(line_no, len(ec)))
            if reference_id not in reference_ids:
                reference_ids.append(reference_id)

            ec_key = ','.join(sorted(reference_ids))

            try:
                ec[ec_key] += 1
            except KeyError:
                ec[ec_key] = 1

            utils.delete_file(temp_file)

    except Exception as e:
        LOG.error("Error: {}".format(str(e)))


    LOG.debug("Process ID: {}, # Unique Reads: {:,}".format(cp.process_id, len(unique_reads)))
    LOG.debug("Process ID: {}, # Reads/Target Duplications: {:,}".format(cp.process_id, same_read_target_counter))
    LOG.debug("Process ID: {}, # Equivalence Classes: {:,}".format(cp.process_id, len(ec)))

    ret = ConvertResults()
    ret.valid_alignments = valid_alignments
    ret.all_alignments = all_alignments
    ret.ec = ec
    ret.unique_reads = unique_reads
    ret.tid_ranges = ranges

    return ret


def process_range_bam(rp):
    """
    NOTE: NOT SURE THIS IS STILL WORKING
    :return:
    """
    LOG.debug('Process ID: {}, Input File: {}'.format(rp.process_id, rp.input_file))

    validate_bam(rp.input_file)

    main_targets = OrderedDict()


    # reference_id: the reference sequence number as defined in the header
    # reference_name: name (None if no AlignmentFile is associated)

    # all the haplotypes
    haplotypes = set()

    # a lookup of reference_ids to main_targets (Ensembl IDs)
    reference_id_to_main_target = {}

    ranges = {}

    all_alignments = 0
    valid_alignments = 0
    reference_id = None
    #temp_name = os.path.join(rp.temp_dir, '_bam2ec.')

    try:
        tmp = {}
        alignment_file = pysam.AlignmentFile(rp.input_file, "rb")
        for target in alignment_file.references:

            idx_underscore = target.rfind('_')
            if idx_underscore > 0:
                main_target = target[:idx_underscore]
            else:
                main_target = target

            if main_target not in tmp:
                tmp[main_target] = main_target
        main_targets_tmp = sorted(tmp.keys())
        for t in main_targets_tmp:
            main_targets[t] = len(main_targets)

        while True:
            alignment = alignment_file.__next__()

            all_alignments += 1

            # if alignment.flag == 4 or alignment.is_unmapped:
            if alignment.is_unmapped:
                continue

            # added for paired end files
            if alignment.is_paired:
                if alignment.is_read2 or not alignment.is_proper_pair or alignment.reference_id != alignment.next_reference_id or alignment.next_reference_start < 0:
                    continue

            valid_alignments += 1

            # reference_sequence_name = Column 3 from file, the Reference NAME (EnsemblID_Haplotype)
            # reference_id = the reference_id, which is 0 or a positive integer mapping to entries
            #       within the sequence dictionary in the header section of a BAM file
            # main_target = the Ensembl id of the transcript

            reference_sequence_name = alignment.reference_name
            reference_id = str(alignment.reference_id)

            idx_underscore = reference_sequence_name.rfind('_')

            if idx_underscore > 0:
                main_target = reference_sequence_name[:idx_underscore]

                if rp.track_ranges:
                    min_max = ranges.get(reference_id, (100000000000, -1))
                    n = min(min_max[0], alignment.reference_start)
                    x = max(min_max[1], alignment.reference_start)
                    ranges[reference_id] = (n,x)

                reference_id_to_main_target[reference_id] = main_target

                try:
                    haplotype = reference_sequence_name[idx_underscore+1:]
                    haplotypes.add(haplotype)
                except:
                    LOG.error('Unable to parse Haplotype from {}'.format(reference_sequence_name))
                    return
            else:
                main_target = reference_sequence_name
                if rp.track_ranges:
                    min_max = ranges.get(reference_id, (100000000000, -1))
                    n = min(min_max[0], alignment.reference_start)
                    x = max(min_max[1], alignment.reference_start)
                    ranges[reference_id] = (n, x)

                reference_id_to_main_target[reference_id] = main_target

                haplotypes.add('')

            if all_alignments % 100000 == 0:
                LOG.debug("Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}".format(rp.process_id, rp.input_file, valid_alignments, all_alignments))

    except StopIteration:
        LOG.info(
            "DONE Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}".format(
                rp.process_id, rp.input_file, valid_alignments, all_alignments))
    except Exception as e:
        print('WHOA'*100)
        print(e)

    haplotypes = sorted(list(haplotypes))

    LOG.debug("# Main Targets: {:,}".format(len(main_targets)))
    LOG.debug("# Haplotypes: {:,}".format(len(haplotypes)))

    ret = RangeResults()
    ret.main_targets = main_targets
    ret.haplotypes = haplotypes
    ret.tid_ranges = ranges

    return ret


def wrapper_convert(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    #print str(args)
    return process_convert_bam(*args)


def wrapper_range(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    #print(str(args))
    return process_range_bam(*args)


def convert(bam_filename, ec_filename, emase_filename, num_chunks=0, number_processes=-1, temp_dir=None, range_filename=None, sample=None, target_filename=None):
    """
    """
    LOG.debug('Parameters')
    LOG.debug('-------------------------------------------')
    LOG.debug('BAM file: {}'.format(bam_filename))
    LOG.debug('EC file: {}'.format(ec_filename))
    LOG.debug('EMASE file: {}'.format(emase_filename))
    LOG.debug('chunks: {}'.format(num_chunks))
    LOG.debug('processes: {}'.format(number_processes))
    LOG.debug('Temp directory: {}'.format(temp_dir))
    LOG.debug('Range file: {}'.format(range_filename))
    LOG.debug('Target file: {}'.format(target_filename))
    LOG.debug('Sample: {}'.format(sample))
    LOG.debug('-------------------------------------------')

    start_time = time.time()

    if number_processes <= 0:
        num_processes = multiprocessing.cpu_count()
    else:
        num_processes = number_processes

    if num_chunks <= 0:
        num_chunks = num_processes
    else:
        if num_chunks > 1000:
            LOG.info("Modifying number of chunks from {} to 1000".format(num_chunks))
            num_chunks = 1000

    # get a sane number of processes
    if num_chunks < num_processes:
        num_processes = num_chunks

    if not temp_dir:
        if emase_filename:
            temp_dir = os.path.dirname(emase_filename)
        if ec_filename:
            temp_dir = os.path.dirname(ec_filename)

    if sample is None:
        sample = os.path.basename(bam_filename)
        LOG.info("Sample not supplied, using filename: {}".format(sample))
    else:
        sample = sample.encode('ascii', 'ignore')

    LOG.info("Parsing file information ...")
    temp_time = time.time()

    alignment_file = pysam.AlignmentFile(bam_filename)
    main_targets = OrderedDict()
    main_targets_list = []
    all_targets_list = []
    target_idx_to_main_target = {}
    haplotypes = set()

    # create main targets from target file
    # target file will be a list of main targets and should always be greater than the number
    # of main targets in bam file
    if target_filename:
        main_targets = utils.parse_targets(target_filename)

        if len(main_targets) == 0:
            LOG.error("Unable to parse target file")
            sys.exit(-1)

        for (k, v) in iteritems(main_targets):
            main_targets_list.append(k)

    #
    for idx, reference_sequence_name in enumerate(alignment_file.references):
        tid = str(alignment_file.get_tid(reference_sequence_name))
        idx_underscore = reference_sequence_name.rfind('_')

        if idx_underscore > 0:
            target = reference_sequence_name[:idx_underscore]
            haplotype = reference_sequence_name[idx_underscore + 1:]
        else:
            target = reference_sequence_name
            haplotype = ''

        target_idx_to_main_target[tid] = target
        haplotypes.add(haplotype)

        if target not in main_targets:
            main_targets[target] = len(main_targets)
            main_targets_list.append(target)

        all_targets_list.append(reference_sequence_name)

    haplotypes = sorted(list(haplotypes))
    haplotypes_idx = {h: idx for idx, h in enumerate(haplotypes)}

    main_target_lengths = np.zeros((len(main_targets), len(haplotypes)), dtype=np.int32)

    alignment_file.close()
    alignment_file = pysam.AlignmentFile(bam_filename)

    #
    # Due to pysam limitation (at least on 0.10.0) we have to LOOP THROUGH TWICE
    # because referencing both alignment_file.lengths and
    # because referencing both alignment_file.references
    #
    for idx, length in enumerate(alignment_file.lengths):
        reference_sequence_name = all_targets_list[idx]
        #print('reference_sequence_name=', reference_sequence_name)

        idx_underscore = reference_sequence_name.rfind('_')

        if idx_underscore > 0:
            target = reference_sequence_name[:idx_underscore]
            haplotype = reference_sequence_name[idx_underscore + 1:]
        else:
            target = reference_sequence_name
            haplotype = ''

        #print('target=', target)
        #print('main_targets[target]=', main_targets[target])
        #print('haplotype=', haplotype)
        #print('haplotypes_idx[haplotype]=', haplotypes_idx[haplotype])

        main_target_lengths[main_targets[target], haplotypes_idx[haplotype]] = length

    LOG.info("File parsed in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                        utils.format_time(start_time, time.time())))

    all_params = []

    LOG.info("Calculating {:,} chunks".format(num_chunks))
    temp_time = time.time()
    chunks = calculate_chunks(bam_filename, num_chunks)
    LOG.info("{:,} chunks calculated in {}, total time: {}".format(len(chunks),
                                                                   utils.format_time(temp_time, time.time()),
                                                                   utils.format_time(start_time, time.time())))
    pid = 0
    for temp_chunk_ids in utils.partition([idx for idx in xrange(num_chunks)], num_processes):
        params = ConvertParams()
        params.input_file = bam_filename
        params.temp_dir = temp_dir
        params.track_ranges = (range_filename is not None)

        for x, cid in enumerate(temp_chunk_ids):
            params.process_id = pid
            params.data.append((cid, chunks[cid]))

        pid += 1
        all_params.append(params)
        LOG.debug('params = {}'.format(str(params)))

    final = ConvertResults()
    final.valid_alignments = 0
    final.all_alignments = 0
    final.ec = OrderedDict()
    final.haplotypes = set()
    final.target_idx_to_main_target = {}
    final.unique_reads = {}
    final.tid_ranges = {}


    temp_time = time.time()
    LOG.info("Starting {} processes ...".format(num_processes))

    args = zip(all_params)
    pool = multiprocessing.Pool(num_processes)

    # KEY = ec_key, VALUE = inserted order
    ec_idx = {}

    for idx, result in enumerate(pool.imap(wrapper_convert, args)):
        LOG.info("Process {} done out of {}, combining result".format(idx + 1, len(all_params)))
        if not final.init:
            final = result
            final.init = True

            # k = ec value, v = count
            for (k, v) in iteritems(final.ec):
                ec_idx[k] = len(ec_idx)

        else:
            # combine ec
            LOG.debug("CHUNK {}: # Result Equivalence Classes: {:,}".format(idx, len(result.ec)))
            for (k, v) in iteritems(result.ec):
                if k in final.ec:
                    final.ec[k] += v
                else:
                    final.ec[k] = v
                    ec_idx[k] = len(ec_idx)
            LOG.debug("CHUNK {}: # Total Equivalence Classes: {:,}".format(idx, len(final.ec)))

            # unique reads
            LOG.debug("CHUNK {}: # Result Unique Reads: {:,}".format(idx, len(result.unique_reads)))
            for (k, v) in iteritems(result.unique_reads):
                if k in final.unique_reads:
                    final.unique_reads[k] += v
                else:
                    final.unique_reads[k] = v
            LOG.debug("CHUNK {}: # Total Unique Reads: {:,}".format(idx, len(final.unique_reads)))

            final.valid_alignments += result.valid_alignments
            final.all_alignments += result.all_alignments

            if range_filename:
                # tid_stats
                for (k, v) in iteritems(result.tid_ranges):
                    if k in final.tid_ranges:
                        n = final.tid_ranges[k][0]
                        x = final.tid_ranges[k][1]
                        final.tid_ranges[k] = (min(v[0], n), max(v[1], x))
                    else:
                        final.tid_ranges[k] = v

        LOG.debug("CHUNK {}: results combined in {}, total time: {}".format(idx, utils.format_time(temp_time, time.time()),
                 utils.format_time(start_time, time.time())))

    LOG.info("All results combined in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
             utils.format_time(start_time, time.time())))

    LOG.info("# Valid Alignments: {:,}".format(final.valid_alignments))
    LOG.info("# Main Targets: {:,}".format(len(main_targets)))
    LOG.info("# Haplotypes: {:,}".format(len(haplotypes)))
    LOG.info("# Equivalence Classes: {:,}".format(len(final.ec)))
    LOG.info("# Unique Reads: {:,}".format(len(final.unique_reads)))

    if range_filename:
        # tid_stats
        with open(range_filename, "w") as fw:
            fw.write("#\t")
            fw.write("\t".join(haplotypes))
            fw.write("\n")

            for main_target in main_targets:
                fw.write(main_target)
                fw.write("\t")

                vals = []

                for haplotype in haplotypes:
                    if len(haplotype) == 0:
                        read_transcript = main_target
                    else:
                        read_transcript = '{}_{}'.format(main_target, haplotype)

                    read_transcript_idx = str(alignment_file.gettid(read_transcript))

                    try:
                        min_max = final.tid_ranges[read_transcript_idx]
                        if min_max[0] == 100000000000 and min_max[1] == -1:
                            vals.append('0')
                        else:
                            vals.append(str(min_max[1] - min_max[0] + 1))
                    except KeyError as ke:
                        vals.append('0')

                fw.write("\t".join(vals))
                fw.write("\n")

    try:
        temp_time = time.time()
        LOG.info('Constructing temp APM structure...')

        # (transcripts, haplotypes, ec)
        new_shape = (len(main_targets),
                     len(haplotypes),
                     len(final.ec))

        LOG.debug('Shape={}'.format(new_shape))

        ec_ids = [x for x in xrange(0, len(final.ec))]
        ec_arr = [[] for _ in xrange(0, len(haplotypes))]
        target_arr = [[] for _ in xrange(0, len(haplotypes))]

        data = []

        # k = comma seperated string of tids
        # v = the count
        for (k, v) in iteritems(final.ec):
            arr_target_idx = k.split(",")

            # get the main targets by name
            temp_main_targets = set()
            for idx in arr_target_idx:
                temp_main_targets.add(target_idx_to_main_target[idx])

            # loop through the targets and haplotypes to get the bits
            for main_target in temp_main_targets:
                # main_target is not an index, but a value like 'ENMUST..001'

                for i, hap in enumerate(haplotypes):
                    if len(hap) == 0:
                        # leaving as 'ENMUST..001'
                        read_transcript = main_target
                    else:
                        # making 'ENMUST..001_A'
                        read_transcript = '{}_{}'.format(main_target, hap)

                    # get the numerical tid corresponding to read_transcript
                    read_transcript_idx = str(alignment_file.gettid(read_transcript))

                    if read_transcript_idx in arr_target_idx:
                        #LOG.debug("{}\t{}\t{}".format(final.ec_idx[k], final.main_targets[main_target], i))

                        # main_targets[main_target] = idx of main target
                        # i = the haplotype
                        # ec_idx[k] = index of ec

                        ec_arr[i].append(ec_idx[k])
                        target_arr[i].append(main_targets[main_target])

                        #print('i = ', i)
                        #print('ec_arr[i] appending ', ec_idx[k])
                        #print('target_arr[i] appending ', main_targets[main_target])

            data.append(v)

        apm = APM(shape=new_shape,
                  haplotype_names=haplotypes,
                  locus_names=main_targets.keys(),
                  read_names=ec_ids,
                  sample_names=[sample])

        for h in xrange(0, len(haplotypes)):
            d = np.ones(len(ec_arr[h]))
            apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])),
                                     shape=(len(final.ec),
                                            len(main_targets)))

        LOG.debug('Constructing CRS...')
        LOG.debug('CRS dimensions: {:,} x {:,}'.format(len(final.ec), 1))

        #LOG.info('len data={}'.format(len(data)))
        #LOG.info('data={}'.format(data[-10:]))

        apm.count = csc_matrix(np.matrix(data).T)

        LOG.info("APM Created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                            utils.format_time(start_time, time.time())))

        if emase_filename:
            try:
                os.remove(emase_filename)
            except OSError:
                pass

            temp_time = time.time()
            LOG.info("Flushing to disk...")
            apm.finalize()
            apm.save(emase_filename, title='bam2ec')
            LOG.info("{} created in {}, total time: {}".format(emase_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))
        if ec_filename:
            LOG.debug("Creating summary matrix...")

            try:
                os.remove(ec_filename)
            except OSError:
                pass

            temp_time = time.time()
            num_haps = len(haplotypes)
            summat = apm.data[0]
            for h in xrange(1, num_haps):
                summat = summat + ((2 ** h) * apm.data[h])

            LOG.debug('summat.sum = {}'.format(summat.sum()))
            LOG.debug('summat.max = {}'.format(summat.max()))
            LOG.debug('summat = {}'.format(summat))

            LOG.info("Matrix created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                   utils.format_time(start_time, time.time())))

            temp_time = time.time()
            LOG.info("Generating BIN file...")

            with open(ec_filename, 'wb') as f:
                # format
                f.write(pack('<i', 2))
                LOG.info("FORMAT: 2")

                #
                # SECTION: HAPLOTYPES
                #     [# of HAPLOTYPES = H]
                #     [length of HAPLOTYPE 1 text][HAPLOTYPE 1 text]
                #     ...
                #     [length of HAPLOTYPE H text][HAPLOTYPE H text]
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

                LOG.info("NUMBER OF HAPLOTYPES: {:,}".format(len(haplotypes)))
                f.write(pack('<i', len(haplotypes)))
                for idx, hap in enumerate(haplotypes):
                    #LOG.debug("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
                    f.write(pack('<i', len(hap)))
                    f.write(pack('<{}s'.format(len(hap)), hap))

                #
                # SECTION: TARGETS
                #     [# of TARGETS = T]
                #     [length TARGET 1 text][TARGET 1 text][HAP 1 length] ... [HAP H length]
                #     ...
                #     [length TARGET T text][TARGET T text][HAP 1 length] ... [HAP H length]
                #
                # Example:
                #     80000
                #     18 ENSMUST00000156068 234
                #     18 ENSMUST00000209341 1054
                #     ...
                #     18 ENSMUST00000778019 1900
                #

                LOG.info("NUMBER OF TARGETS: {:,}".format(len(main_targets)))
                f.write(pack('<i', len(main_targets)))
                for (main_target, idx) in iteritems(main_targets):
                    f.write(pack('<i', len(main_target)))
                    f.write(pack('<{}s'.format(len(main_target)), main_target))

                    lengths = []

                    for idx_hap, hap in enumerate(haplotypes):
                        length = main_target_lengths[idx, idx_hap]
                        f.write(pack('<i', length))
                        #lengths.append(str(length))

                    #LOG.debug("#{:,} --> {:,}\t{}\t{}\t".format(idx, len(main_target), main_target, '\t'.join(lengths)))

                #
                # SECTION: CRS
                #     [# of CRS = 1]
                #     [length of CR 1 text][CR 1 text]
                #
                # Example:
                #     1
                #     8 SAMPLEID
                #

                temp_sample = sample.encode('utf-8')

                LOG.info("FILTERED CRS: 1")
                f.write(pack('<i', 1))
                f.write(pack('<i', len(temp_sample)))
                f.write(pack('<{}s'.format(len(temp_sample)), temp_sample))


                #
                # SECTION: ALIGNMENT MAPPINGS ("A" Matrix)
                #     [# of ALIGNMENT MAPPINGS (AM) = A]
                #     [EC INDEX][TRANSCRIPT INDEX][HAPLOTYPE flag] (for AM 1)
                #     [EC INDEX][TRANSCRIPT INDEX][HAPLOTYPE flag] (for AM 2)
                #     ...
                #     [EC INDEX][TRANSCRIPT INDEX][HAPLOTYPE flag] (for AM A)
                #
                # NOTE:
                #     HAPLOTYPE flag is an integer that denotes which haplotype
                #     (allele) a read aligns to given an EC. For example, 00, 01,
                #     10, and 11 can specify whether a read aligns to the 1st
                #     and/or 2nd haplotype of a transcript.  These binary numbers
                #     are converted to integers - 0, 1, 2, 3 - and stored as the
                #     haplotype flag.
                #
                # Example:
                #     5000
                #     1 2 4
                #     8 2 1
                #     ...
                #     100 200 8
                #

                LOG.info("Determining mappings...")

                num_mappings = summat.nnz
                summat = coo_matrix(summat)
                summat = summat.tocsr()

                LOG.info("A MATRIX: INDPTR LENGTH {:,}".format(len(summat.indptr)))
                f.write(pack('<i', len(summat.indptr)))

                # NON ZEROS
                LOG.info("A MATRIX: NUMBER OF NON ZERO: {:,}".format(num_mappings))
                f.write(pack('<i', num_mappings))

                # ROW OFFSETS
                LOG.info("A MATRIX: LENGTH INDPTR: {:,}".format(len(summat.indptr)))
                f.write(pack('<{}i'.format(len(summat.indptr)), *summat.indptr))
                # LOG.error(summat.indptr)

                # COLUMNS
                LOG.info("A MATRIX: LENGTH INDICES: {:,}".format(len(summat.indices)))
                f.write(pack('<{}i'.format(len(summat.indices)), *summat.indices))
                # LOG.error(summat.indices)

                # DATA
                LOG.info("A MATRIX: LENGTH DATA: {:,}".format(len(summat.data)))
                f.write(pack('<{}i'.format(len(summat.data)), *summat.data))
                # LOG.error(summat.data)

                #
                # SECTION: "N" Matrix
                #
                # "N" Matrix format is EC (rows) by CRS (columns) with
                # each value being the EC count.
                #
                # Instead of storing a "dense" matrix, we store a "sparse"
                # matrix utilizing Compressed Sparse Column (CSC) format.
                #

                LOG.info("N MATRIX: NUMBER OF EQUIVALENCE CLASSES: {:,}".format(len(final.ec)))
                LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<i', len(apm.count.indptr)))

                # NON ZEROS
                LOG.info("N MATRIX: NUMBER OF NON ZERO: {:,}".format(apm.count.nnz))
                f.write(pack('<i', apm.count.nnz))

                # ROW OFFSETS
                LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<{}i'.format(len(apm.count.indptr)), *apm.count.indptr))
                # LOG.error(apm.count.indptr)

                # COLUMNS
                LOG.info("N MATRIX: LENGTH INDICES: {:,}".format(len(apm.count.indices)))
                f.write(pack('<{}i'.format(len(apm.count.indices)), *apm.count.indices))
                # LOG.error(apm.count.indices)

                # DATA
                LOG.info("N MATRIX: LENGTH DATA: {:,}".format(len(apm.count.data)))
                f.write(pack('<{}i'.format(len(apm.count.data)), *apm.count.data))
                # LOG.error(apm.count.data)


            LOG.info("{} created in {}, total time: {}".format(ec_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

    except KeyboardInterrupt as e:
        LOG.fatal("ERROR: {}".format(str(e)))
        raise Exception(e)


def split_bam(filename, n, directory=None):
    """
    Split a BAM file into ``n`` smaller bam files.

    :param str filename: the name of the BAM file
    :param int n: number of files to chunk into
    :param str directory: output directory, defaults to ``bam_filename`` directory

    :return: names of the files
    """
    start_time = time.time()

    LOG.debug("BAM File: {}".format(filename))
    LOG.debug("Number of Files: {}".format(n))

    if not directory:
        directory = os.path.dirname(filename)

    LOG.debug("Output Directory: {}".format(directory))

    bam_basename = os.path.basename(filename)
    bam_prefixname, bam_extension = os.path.splitext(bam_basename)
    bam_output_temp = os.path.join(directory, bam_prefixname)

    LOG.info("Calculating {:,} chunks...".format(n))
    temp_time = time.time()
    chunks = calculate_chunks(filename, n)
    LOG.info("{:,} chunks calculated in {}, total time: {}".format(len(chunks),
                                                                   utils.format_time(temp_time, time.time()),
                                                                   utils.format_time(start_time, time.time())))

    file_names = []

    for idx, chunk in enumerate(chunks):
        # must create the file
        LOG.debug(chunk)
        new_file = "{}_{}{}".format(bam_output_temp, idx, bam_extension)
        file_names.append(new_file)
        LOG.debug("Creating alignment file: {}".format(new_file))
        chunk_bam_file(filename, new_file, chunk)

    LOG.info("{:,} files created in {}, total time: {}".format(len(chunks),
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

    return file_names


def bytes_from_file_bam(read_filename, write_filename, offset=0, bytes_size=-1):
    """
    Read bytes from a file and append them onto another file.

    :param read_filename: the name of the file to read from
    :param write_filename: the name of the file to write to
    :param offset: the number of bytes to offset (seek)
    :param bytes_size: the number of bytes to read, -1 = to end of file
    :return:
    """
    try:

        with open(read_filename, "rb") as fr:
            if offset > 0:
                fr.seek(offset)

            if bytes_size > 0:
                data = fr.read(bytes_size)
            else:
                data = fr.read()

            mode = 'r+b'
            size = os.path.getsize(write_filename)

            with open(write_filename, mode) as fw:
                fw.seek(size - 28)
                temp = fw.read()
                if temp == BAM_EOF:
                    fw.seek(size - 28)
                # else:
                #    fw.seek(size)
                fw.write(data)
    except Exception as e:
        LOG.error('ERROR: {}'.format(str(e)))



"""
import pysam
from Bio import bgzf
fhin = "../bam2ec/data/smaller_file/bowtie.bam"
fhout = open("tells", "w")
try:
    aln_file = pysam.AlignmentFile(fhin)


    i = 1
    while True:
        tell = aln_file.tell()
        aln = aln_file.__next__()
        data = bgzf.split_virtual_offset(tell)
        fhout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i, data[0], data[1], tell, aln.qname, aln.tid))
        i += 1
except StopIteration:
    pass
fhout.close()
"""


def calculate_chunks(filename, num_chunks):
    """
    Calculate the boundaries in the BAM file and partition into chunks.

    :param str filename: name of the BAM file
    :param int num_chunks: number of chunks to partition the boundaries into
    :return: a list of tuples containing the start and end boundaries
    """
    if num_chunks <= 0:
        raise ValueError("The number of chunks to calculate should be >= 1")

    header_size = get_header_size(filename)

    if num_chunks == 1:
        #aln_file = pysam.AlignmentFile(filename)
        #header_size = bgzf.split_virtual_offset(aln_file.tell())[0]
        #aln_file.close()

        pr = ParseRecord(header_size, 0, 0, header_size, -1, 0, 0)
        return [pr]

    try:
        f = open(filename, 'rb')
        # get all the block start offsets
        block_offsets = []
        decompressed_lengths = []
        i = 0

        for values in FastBgzfBlocks(f):
        #for values in bgzf.BgzfBlocks(f):
            # make sure that the first block is greater than header size
            # we don't want to split mid header even though this should only happen
            # on large values for num_chunks
            if header_size > values[0]:
                continue

            block_offsets.append(values[0])
            decompressed_lengths.append(values[3])

            if i % 50000 == 0:
                LOG.debug('Block {}'.format(i))
            i += 1

        # partition the starts into manageable chunks
        div, mod = divmod(len(block_offsets), num_chunks)

        aln_file = pysam.AlignmentFile(filename)
        header_size = bgzf.split_virtual_offset(aln_file.tell())[0]
        partitioned_offsets = [(header_size, 0)]

        for i in xrange(1, num_chunks):
            index = div * i + min(i, mod)
            virtual_offset = bgzf.make_virtual_offset(block_offsets[index], 0)
            aln_file.seek(virtual_offset)
            aln = aln_file.__next__()
            aln_qname = aln.query_name

            # added because some query_name had spaces
            i = aln_qname.find(' ')
            if i > 0:
                aln_qname = aln_qname[:i]

            aln_first_qname = aln_qname

            while aln_first_qname == aln_qname:

                virtual_offset = aln_file.tell()
                aln = aln_file.__next__()
                aln_qname = aln.query_name
                i = aln_qname.find(' ')
                if i > 0:
                    aln_qname = aln_qname[:i]

            partitioned_offsets.append(bgzf.split_virtual_offset(virtual_offset))

        aln_file.close()

        # now let's calculate beginning and ends
        params = []

        for i, offset in enumerate(partitioned_offsets):
            #print '{} => {}'.format(i, offset)

            index = block_offsets.index(partitioned_offsets[i][0])
            begin_read_offset = 0
            begin_read_size = 0
            file_offset = 0
            file_bytes = 0
            end_read_offset = 0
            end_read_size = 0

            if i == 0:
                # first
                begin_read_offset = 0
                begin_read_size = 0
                file_offset = block_offsets[index]
                #print 'file_offset=', file_offset
                file_bytes = partitioned_offsets[i + 1][0] - file_offset
                #print 'file_bytes=', file_bytes
                end_read_offset = bgzf.make_virtual_offset(partitioned_offsets[i + 1][0], 0)
                end_read_size = partitioned_offsets[i + 1][1]
            elif i == num_chunks - 1:
                # last
                begin_read_offset = bgzf.make_virtual_offset(partitioned_offsets[i][0], partitioned_offsets[i][1])
                begin_read_size = decompressed_lengths[index] - partitioned_offsets[i][1]
                file_offset = block_offsets[index + 1]
                file_bytes = -1
                end_read_offset = 0
                end_read_size = 0
            else:
                # all others
                if offset[1] == 0:
                    # bgzf boundary
                    LOG.fatal('****************HUH')
                    return

                begin_read_offset = bgzf.make_virtual_offset(partitioned_offsets[i][0], partitioned_offsets[i][1])
                begin_read_size = decompressed_lengths[index] - partitioned_offsets[i][1]
                file_offset = block_offsets[index + 1]
                file_bytes = partitioned_offsets[i + 1][0] - file_offset
                end_read_offset = bgzf.make_virtual_offset(partitioned_offsets[i + 1][0], 0)
                end_read_size = partitioned_offsets[i + 1][1]

            pr = ParseRecord(header_size, begin_read_offset, begin_read_size, file_offset, file_bytes, end_read_offset,
                             end_read_size)
            params.append(pr)

        return params

    except Exception as e:
        LOG.debug('calculate_chunks error: {}'.format(str(e)))



def truncate_bam_file(fname):
    """
    Remove the EOF from BGZF/BAM file.

    Does not check if the EOF is present or not.

    :param fname: the name of the BZF/BAM file
    :return:
    """
    utils.truncate_file(fname, 28)


def FastBgzfBlocks(handle):
    """
    Faster version of bgzf.BgzfBlocks

    :param handle: the handle to the BGZF file
    :return: tuple of start offset, block length, data offset, data length
    """
    data_start = 0
    while True:
        start_offset = handle.tell()
        # This may raise StopIteration which is perfect here
        block_length, data_len = _quick_bgzf_load(handle)
        yield start_offset, block_length, data_start, data_len
        data_start += data_len


def _quick_bgzf_load(handle):
    """
    Quicker version of bgzf._bgzf_load.  No decompressing of BGZF data.  Just getting meta information.
    """
    magic = handle.read(4)

    if not magic:
        raise StopIteration

    if magic != bgzf._bgzf_magic:
        raise ValueError("A BGZF block should start with %r, not %r; handle.tell() now says %r" % (bgzf._bgzf_magic, magic, handle.tell()))

    gzip_mod_time, gzip_extra_flags, gzip_os, extra_len = struct.unpack("<LBBH", handle.read(8))
    block_size = None
    x_len = 0

    while x_len < extra_len:
        subfield_id = handle.read(2)
        subfield_len = struct.unpack("<H", handle.read(2))[0]
        subfield_data = handle.read(subfield_len)
        x_len += subfield_len + 4
        if subfield_id == bgzf._bytes_BC:
            assert subfield_len == 2, "Wrong BC payload length"
            assert block_size is None, "Two BC subfields?"
            block_size = struct.unpack("<H", subfield_data)[0] + 1
    assert x_len == extra_len, (x_len, extra_len)
    assert block_size is not None, "Missing BC, this isn't a BGZF file!"
    deflate_size = block_size - 1 - extra_len - 19
    handle.seek(handle.tell() + deflate_size)
    expected_crc = handle.read(4)
    expected_size = struct.unpack("<I", handle.read(4))[0]
    return block_size, expected_size


def generate_bam_ranges(input_files, range_filename, target_filename=None, temp_dir=None):
    """
    we might want to check if the number of files is less than the number of cpus

    if that's the case, we can split the input files and process in threads

    :param input_files:
    :param temp_dirs:
    :return:
    """
    start_time = time.time()

    num_processes = multiprocessing.cpu_count()

    """
    if len(input_files) < num_processes:
        # split files, we could probably do better but...

        files = []
        for filename in input_files:
            size = os.path.getsize(filename)
            files.append((size, filename))

        while len(files) < num_processes:
            print('len(files) < num_processes')
            print('{} < {}'.format(len(files), num_processes))
            # sort files by size descending
            input_files = sorted(files, key=lambda x: x[0], reverse=True)
            print(' splitting')
            new_names = split_bam(input_files[0][1], 2)

            files = input_files[1:]
            files.append((os.path.getsize(new_names[0]), new_names[0]))
            files.append((os.path.getsize(new_names[1]), new_names[1]))

        input_files = files

    for f in input_files:
        print(f)

    """

    if not temp_dir:
        temp_dir = os.getcwd()

    temp_files = []

    # sort files based upon size, take the largest and split it
    files = []
    for filename in input_files:
        size = os.path.getsize(filename)
        files.append((size, filename))

    if len(input_files) < num_processes:
        LOG.info("Preparing files...")
        temp_time = time.time()

        input_files = sorted(files, key=lambda x: x[0], reverse=True)
        new_names = split_bam(input_files[0][1], num_processes - len(files) + 1, directory=temp_dir)

        files = input_files[1:]
        for f in new_names:
            files.append((os.path.getsize(f), f))
            temp_files.append(f)

        LOG.info("{:,} files prepared in {}, total time: {}".format(len(input_files),
                                                                   utils.format_time(temp_time, time.time()),
                                                                   utils.format_time(start_time, time.time())))
    input_files = files

    all_params = []
    for i, f in enumerate(input_files):
        params = RangeParams()
        params.process_id = i
        params.input_file = f[1]
        print (params.input_file)

        params.target_file = target_filename
        params.temp_dir = temp_dir
        all_params.append(params)

    LOG.info("Starting {} processes ...".format(num_processes))

    temp_time = time.time()
    args = zip(all_params)
    pool = multiprocessing.Pool(num_processes)
    results = pool.map(wrapper_range, args)

    LOG.info("All processes done in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

    LOG.info("Combining {} results ...".format(len(results)))
    temp_time = time.time()

    final = ConvertResults()
    final.haplotypes = set()
    final.main_targets = OrderedDict()
    final.tid_ranges = {}

    alignment_file = pysam.AlignmentFile(input_files[0][1])

    # parse results
    for idx, result in enumerate(results):
        if not final.init:
            final = result
            final.init = True
        else:
            # assuming haplotypes and main_targets are same
            if final.main_targets != result.main_targets:
                raise ValueError('Error: main targets do not match')

            if final.haplotypes != result.haplotypes:
                raise ValueError('Error: haplotypes do not match')

            # tid_stats
            for (k, v) in iteritems(result.tid_ranges):
                if k in final.tid_ranges:
                    n = final.tid_ranges[k][0]
                    x = final.tid_ranges[k][1]
                    final.tid_ranges[k] = (min(v[0], n), max(v[1], x))
                else:
                    final.tid_ranges[k] = v

        LOG.debug("CHUNK {}: results combined in {}, total time: {}".format(idx, utils.format_time(temp_time, time.time()),
                 utils.format_time(start_time, time.time())))

    LOG.info("All results combined in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
             utils.format_time(start_time, time.time())))

    LOG.info("# Main Targets: {:,}".format(len(final.main_targets)))
    LOG.info("# Haplotypes: {:,}".format(len(final.haplotypes)))

    if temp_files:
        LOG.info("Cleaning up temporary files...")

        for t in temp_files:
            utils.delete_file(t)

    LOG.info("Generating range file {} ...".format(range_filename))

    temp_time = time.time()

    with open(range_filename, "w") as fw:
        fw.write("#\t")
        fw.write("\t".join(final.haplotypes))
        fw.write("\n")

        for main_target in final.main_targets:
            fw.write(main_target)
            fw.write("\t")

            vals = []

            for haplotype in final.haplotypes:
                read_transcript = '{}_{}'.format(main_target, haplotype)
                read_transcript_idx = str(alignment_file.gettid(read_transcript))

                try:
                    min_max = final.tid_ranges[read_transcript_idx]
                    if min_max[0] == 100000000000 and min_max[1] == -1:
                        vals.append('0')
                    else:
                        vals.append(str(min_max[1] - min_max[0] + 1))
                except KeyError as ke:
                    vals.append('0')

            fw.write("\t".join(vals))
            fw.write("\n")

    LOG.info("{} created in {}, total time: {}".format(range_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))
