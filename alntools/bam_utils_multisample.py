# -*- coding: utf-8 -*-
from collections import OrderedDict, namedtuple, defaultdict
from struct import pack

import glob
import multiprocessing
import os
import struct
import sys
import time

from Bio import bgzf
from scipy.sparse import coo_matrix

import numpy as np
import pysam

from . import utils
from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

try:
    xrange
except NameError:
    xrange = range


LOG = utils.get_logger()
BAM_HEADER = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
BAM_EOF = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

parse_fields = ["header_size", "begin_read_offset", "begin_read_size", "file_offset", "file_bytes", "end_read_offset",
                "end_read_size"]
ParseRecord = namedtuple("ParseRecord", parse_fields)


class ConvertParams(object):
    """
    # each core needs
    # - name of alignment file
    # - header size
    # - target file
    # - list of files to create and work on
    #   - idx, vo_start, vo_end

    """
    slots = ['input_file', 'target_file', 'temp_dir', 'process_id', 'data', 'emase', 'track_ranges']

    def __init__(self):
        self.input_file = None
        self.target_file = None
        self.temp_dir = None
        self.process_id = None
        self.data = [] # tuple of (idx, ParseRecord)
        self.emase = False
        self.track_ranges = False

    def __str__(self):
        return "Input: {}\nProcess ID: {}\nData: {}".format(self.input_file, self.process_id, self.data)


class ConvertResults(object):
    slots = ['main_targets', 'valid_alignments', 'all_alignments', 'ec', 'ec_idx', 'haplotypes', 'target_idx_to_main_target', 'unique_reads', 'init', 'tid_ranges', 'CRS', 'CRS_idx']

    def __init__(self):
        self.main_targets = None
        self.valid_alignments = None
        self.all_alignments = None
        self.ec = None
        self.ec_idx = None
        self.haplotypes = None
        self.target_idx_to_main_target = None
        self.unique_reads = None
        self.init = False
        self.tid_ranges = None
        self.CRS = None
        self.CRS_idx = None


class RangeParams(object):
    """
    # each core needs
    # - name of alignment file
    # - header size
    # - target file
    # - list of files to create and work on
    #   - idx, vo_start, vo_end

    """
    slots = ['input_file', 'target_file', 'temp_dir', 'process_id']

    def __init__(self):
        self.input_file = None
        self.target_file = None
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


def ddict2dict(d):
    for k, v in d.items():
        if isinstance(v, dict):
            d[k] = ddict2dict(v)
    return dict(d)


def process_convert_bam(cp):
    """

    :return:
    """
    LOG.debug('Process ID: {}, Input File: {}'.format(cp.process_id, cp.input_file))

    if cp.target_file:
        LOG.debug('Process ID: {}, Target File: {}'.format(cp.process_id, cp.target_file))

    if cp.emase:
        LOG.debug('Process ID: {}, Emase format requested'.format(cp.process_id))

    validate_bam(cp.input_file)

    main_targets = OrderedDict()

    if cp.target_file:
        main_targets = utils.parse_targets(cp.target_file)
        if len(main_targets) == 0:
            LOG.error("Unable to parse target file")
            sys.exit(-1)
    else:
        tmp = {}
        alignment_file = pysam.AlignmentFile(cp.input_file, "rb")
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
        alignment_file.close()


    # reference_id: the reference sequence number as defined in the header
    # reference_name: name (None if no AlignmentFile is associated)

    # ec = equivalence class
    #      the KEY is a comma separated string of reference_ids
    #      the VALUE is the number of times this equivalence class has appeared
    ec = defaultdict(lambda: defaultdict(int))

    # ec_idx = lookup to ec
    #          the KEY is a comma separated string of reference_ids
    #          the VALUE is a number specifying the insertion order of the KEY value in ec
    ec_idx = {}

    # all the haplotypes
    haplotypes = set()

    # a lookup of reference_ids to main_targets (Ensembl IDs)
    reference_id_to_main_target = {}

    # unique reads
    unique_reads = {}

    # times encountering new read id
    read_id_switch_counter = 0

    same_read_target_counter = 0

    ranges = {}

    #pid = os.getpid()

    all_alignments = 0
    valid_alignments = 0
    ec_key = None
    temp_name = os.path.join(cp.temp_dir, '_bam2ec.')
    generic_haplotype = '0'

    reference_id = None
    reference_ids = []

    try:
        alignment_file = pysam.AlignmentFile(cp.input_file)

        # query template name
        query_name = None

        while True:
            alignment = alignment_file.next()

            all_alignments += 1

            # if alignment.flag == 4 or alignment.is_unmapped:
            if alignment.is_unmapped:
                continue

            # added for paired end files
            if alignment.is_paired:
                if alignment.is_read2 or not alignment.is_proper_pair or alignment.reference_id != alignment.next_reference_id or alignment.next_reference_start < 0:
                    continue

            '''
            print('alignment.query_name={}'.format(alignment.query_name))
            print('alignment.reference_id={}'.format(alignment.reference_id))
            print('alignment.reference_name={}'.format(alignment.reference_name))
            print('alignment.is_paired={}'.format(alignment.is_paired))
            print('alignment.is_proper_pair={}'.format(alignment.is_proper_pair))
            print('alignment.next_reference_id={}'.format(alignment.next_reference_id))
            print('alignment.next_reference_start={}'.format(alignment.next_reference_start))
            print('alignment.is_unmapped={}'.format(alignment.is_unmapped))
            '''

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

                if cp.track_ranges:
                    min_max = ranges.get(reference_id, (100000000000, -1))
                    n = min(min_max[0], alignment.reference_start)
                    x = max(min_max[1], alignment.reference_start)
                    ranges[reference_id] = (n,x)

                if cp.target_file:
                    if main_target not in main_targets:
                        LOG.error("Unexpected target found in BAM file: {}".format(main_target))
                        sys.exit(-1)

                reference_id_to_main_target[reference_id] = main_target

                try:
                    haplotype = reference_sequence_name[idx_underscore+1:]
                    haplotypes.add(haplotype)
                except:
                    LOG.error('Unable to parse Haplotype from {}'.format(reference_sequence_name))
                    return
            else:
                main_target = reference_sequence_name
                if cp.track_ranges:
                    min_max = ranges.get(reference_id, (100000000000, -1))
                    n = min(min_max[0], alignment.reference_start)
                    x = max(min_max[1], alignment.reference_start)
                    ranges[reference_id] = (n, x)

                if cp.target_file:
                    if main_target not in main_targets:
                        LOG.error("Unexpected target found in BAM file: {}".format(main_target))
                        sys.exit(-1)

                reference_id_to_main_target[reference_id] = main_target

                haplotype = generic_haplotype
                haplotypes.add(haplotype)

            # query_name = Column 1 from file, the Query template NAME

            if query_name is None:
                query_name = alignment.query_name

            # query_name = Column 1 from file, the Query template NAME
            bam_tags = query_name.split('|||')
            orig_query_name = bam_tags[0]
            bam_tag_CR = bam_tags[2]
            bam_tag_CY = bam_tags[4]
            bam_tag_UR = bam_tags[6]
            bam_tag_UY = bam_tags[8]
            bam_tag_BC = bam_tags[10]
            bam_tag_QT = bam_tags[12]

            # set tags and dump new alignment

            try:
                unique_reads[query_name] += 1
            except KeyError:
                unique_reads[query_name] = 1

            alignment_query_name = alignment.query_name

            if query_name != alignment_query_name:
                ec_key = ','.join(sorted(reference_ids))
                ec[ec_key][bam_tag_CR] += 1

                query_name = alignment.query_name

                reference_ids = [reference_id]
                read_id_switch_counter += 1
            else:
                if reference_id not in reference_ids:
                    reference_ids.append(reference_id)
                else:
                    same_read_target_counter += 1

            if all_alignments % 100000 == 0:
                LOG.debug("Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}, with {:,} equivalence classes".format(cp.process_id, cp.input_file, valid_alignments, all_alignments, len(ec)))

    except StopIteration:
        LOG.info(
            "DONE Process ID: {}, File: {}, {:,} valid alignments processed out of {:,}, with {:,} equivalence classes".format(
                cp.process_id, cp.input_file, valid_alignments, all_alignments, len(ec)))

    haplotypes = sorted(list(haplotypes))

    LOG.debug("Process ID: {}, # Unique Reads: {:,}".format(cp.process_id, len(unique_reads)))
    LOG.debug("Process ID: {}, # Reads/Target Duplications: {:,}".format(cp.process_id, same_read_target_counter))
    LOG.debug("Process ID: {}, # Main Targets: {:,}".format(cp.process_id, len(main_targets)))
    LOG.debug("Process ID: {}, # Haplotypes: {:,}".format(cp.process_id, len(haplotypes)))
    LOG.debug("Process ID: {}, # Equivalence Classes: {:,}".format(cp.process_id, len(ec)))

    ret = ConvertResults()
    ret.main_targets = main_targets
    ret.valid_alignments = valid_alignments
    ret.all_alignments = all_alignments
    ret.ec = ddict2dict(ec)
    ret.ec_idx = OrderedDict()
    ret.haplotypes = haplotypes
    ret.target_idx_to_main_target = reference_id_to_main_target
    ret.unique_reads = unique_reads
    ret.tid_ranges = ranges
    ret.CRS = OrderedDict()

    return ret


def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


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

# barcode_file
# read_id
# sequence
# +
# quality


def convert(bam_filename, output_filename, num_chunks=0, target_filename=None, emase=False, temp_dir=None, range_filename=None, minimum_count=-1, number_processes=-1):
    """
    """
    start_time = time.time()

    if os.path.isfile(bam_filename):
        LOG.error('bam file must be a directory')
        return None

    bam_files = glob.glob(os.path.join(bam_filename, "*.bam"))
    if len(bam_files) == 0:
        LOG.error('No bam files found in directory: {}'.format(bam_filename))
        return None

    if number_processes <= 0:
        num_processes = multiprocessing.cpu_count()
        num_processes = min(num_processes, len(bam_files))
    else:
        num_processes = number_processes

    if not temp_dir:
        temp_dir = os.path.dirname(output_filename)


    # each core needs
    # - name of alignment file
    # - target file
    # - list of files to create and work on
    #   - idx, mp_param

    all_params = []

    pid = 0
    for bam_file in bam_files:
        params = ConvertParams()
        params.input_file = bam_file
        params.target_file = target_filename
        params.temp_dir = temp_dir
        params.track_ranges = (range_filename != None)
        params.process_id = pid
        pid += 1
        all_params.append(params)
        LOG.debug('params = {}'.format(str(params)))

    final = ConvertResults()
    final.valid_alignments = 0
    final.all_alignments = 0
    final.ec = OrderedDict()
    final.ec_idx = OrderedDict()
    final.CRS = {}
    final.haplotypes = set()
    final.main_targets = OrderedDict()
    final.target_idx_to_main_target = {}
    final.unique_reads = {}
    final.tid_ranges = {}

    ec_totals = OrderedDict()
    cr_totals = OrderedDict()

    LOG.info("Starting {} processes ...".format(num_processes))

    temp_time = time.time()
    args = zip(all_params)
    pool = multiprocessing.Pool(num_processes)
    results = pool.map(wrapper_convert, args)

    LOG.info("All processes done in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

    LOG.info("Combining {} results ...".format(len(results)))
    temp_time = time.time()

    alignment_file = pysam.AlignmentFile(bam_files[0])

    counter_test = 0

    # parse results
    for idx, result in enumerate(results):

        LOG.info("Combining result #{}".format(idx))

        if not final.init:
            final = result
            final.init = True
            final.ec_idx = OrderedDict()
            final.CRS_idx = OrderedDict()

            for eckey, crdict in result.ec.iteritems():
                final.ec_idx[eckey] = len(final.ec_idx)
                for crkey, count in crdict.iteritems():
                    final.CRS[crkey] = 1

                    if crkey not in final.CRS_idx:
                        final.CRS_idx[crkey] = len(final.CRS_idx)

                    if crkey in cr_totals:
                        cr_totals[crkey] += count
                    else:
                        cr_totals[crkey] = count

                    if eckey in ec_totals:
                        ec_totals[eckey] += count
                    else:
                        ec_totals[eckey] = count
                        counter_test += 1

        else:
            # combine ec and URS
            LOG.debug("CHUNK {}: # Result Equivalence Classes: {:,}".format(idx, len(result.ec)))

            for eckey, crdict in result.ec.iteritems():
                for crkey, count in crdict.iteritems():
                    final.CRS[crkey] = 1

                    if crkey not in final.CRS_idx:
                        final.CRS_idx[crkey] = len(final.CRS_idx)

                    if crkey in cr_totals:
                        cr_totals[crkey] += count
                    else:
                        cr_totals[crkey] = count

                    if eckey in ec_totals:
                        ec_totals[eckey] += count
                    else:
                        ec_totals[eckey] = count

                    if eckey in final.ec:
                        if crkey in final.ec[eckey]:
                            final.ec[eckey][crkey] += count
                        else:
                            final.ec[eckey][crkey] = count
                            counter_test += 1
                    else:
                        final.ec_idx[eckey] = len(final.ec_idx)
                        final.ec[eckey] = {crkey: count}
                        counter_test += 1

            LOG.debug("CHUNK {}: # Total Equivalence Classes: {:,}".format(idx, len(final.ec)))
            LOG.debug("CHUNK {}: # Total CRs: {:,}".format(idx, len(final.CRS)))

            # combine haplotypes
            LOG.debug("CHUNK {}: # Result Haplotypes: {:,}".format(idx, len(result.haplotypes)))
            s1 = set(final.haplotypes)
            s2 = set(result.haplotypes)
            final.haplotypes = sorted(list(s1.union(s2)))
            LOG.debug("CHUNK {}: # Total Haplotypes: {:,}".format(idx, len(final.haplotypes)))

            # combine target_idx_to_main_target
            LOG.debug("CHUNK {}: # Result Main Targets: {:,}".format(idx, len(result.main_targets)))
            for k, v in result.target_idx_to_main_target.iteritems():
                if k not in final.target_idx_to_main_target:
                    final.target_idx_to_main_target[k] = v
            LOG.debug("CHUNK {}: # Total Main Targets: {:,}".format(idx, len(final.main_targets)))

            # unique reads
            LOG.debug("CHUNK {}: # Result Unique Reads: {:,}".format(idx, len(result.unique_reads)))
            for k, v in result.unique_reads.iteritems():
                if k in final.unique_reads:
                    final.unique_reads[k] += v
                else:
                    final.unique_reads[k] = v
            LOG.debug("CHUNK {}: # Total Unique Reads: {:,}".format(idx, len(final.unique_reads)))

            final.valid_alignments += result.valid_alignments
            final.all_alignments += result.all_alignments

            if range_filename:
                # tid_stats
                for k, v in result.tid_ranges.iteritems():
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

    #LOG.info("# Total Alignments: {:,}".format(final.all_alignments))
    LOG.info("# Valid Alignments: {:,}".format(final.valid_alignments))
    LOG.info("# Main Targets: {:,}".format(len(final.main_targets)))
    LOG.info("# Haplotypes: {:,}".format(len(final.haplotypes)))
    LOG.info("# Equivalence Classes: {:,}".format(len(final.ec)))
    LOG.debug("# Equivalence Classes (ec_idx): {:,}".format(len(final.ec_idx)))
    LOG.debug("# Equivalence Class Max Index: {:,}".format(max(final.ec_idx.values())))
    LOG.info("# Unique Reads: {:,}".format(len(final.unique_reads)))

    # filter everything
    LOG.debug("Minimum Count: {:,}".format(minimum_count))

    if minimum_count > 0:
        LOG.info("FILTERING CRS: {:,}".format(len(final.CRS)))
        # find the new CRS
        final.CRS = OrderedDict()
        final.CRS_idx = OrderedDict()

        for CR, CR_total in cr_totals.iteritems():
            if CR_total >= minimum_count:
                final.CRS[CR] = CR_total

                if CR not in final.CRS_idx:
                    final.CRS_idx[CR] = len(final.CRS_idx)


        # remove invalid CRS from ECs
        new_ecs = OrderedDict()
        new_ec_idx = OrderedDict()
        new_ec_totals = {}
        # loop through ecs
        for eckey, crs in final.ec.iteritems():
            # potential new ec

            ec = {}
            total = 0
            # loop through valid CRS and and if valid
            for crkey, crcount in crs.iteritems():

                if crkey in final.CRS:
                    ec[crkey] = crs[crkey]
                    total += crcount

            # only add to new ecs if there is anything
            if len(ec) > 0:
                new_ecs[eckey] = ec
                new_ec_idx[eckey] = len(new_ec_idx)
                new_ec_totals[eckey] = total

        ec_totals = new_ec_totals
        final.ec = new_ecs
        final.ec_idx = new_ec_idx

    #LOG.info("# Total Alignments: {:,}".format(final.all_alignments))
    LOG.info("# Valid Alignments: {:,}".format(final.valid_alignments))
    LOG.info("# Main Targets: {:,}".format(len(final.main_targets)))
    LOG.info("# Haplotypes: {:,}".format(len(final.haplotypes)))
    LOG.info("# Equivalence Classes: {:,}".format(len(final.ec)))
    LOG.debug("# Equivalence Classes (ec_idx): {:,}".format(len(final.ec_idx)))
    LOG.debug("# Equivalence Class Max Index: {:,}".format(max(final.ec_idx.values())))
    LOG.info("# Unique Reads: {:,}".format(len(final.unique_reads)))

    if range_filename:
        # tid_stats
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

    try:
        os.remove(output_filename)
    except OSError:
        pass

    LOG.info("CRS: {:,}".format(len(final.CRS)))
    ec_arr_max = -1
    target_arr_max = -1

    if emase:
        try:
            temp_time = time.time()
            LOG.info('Constructing APM structure...')

            new_shape = (len(final.main_targets), len(final.haplotypes), len(final.ec))

            ec_ids = [x for x in xrange(0, len(final.ec))]

            # ec.values -> the number of times this equivalence class has appeared

            ec_arr = [[] for _ in xrange(0, len(final.haplotypes))]
            target_arr = [[] for _ in xrange(0, len(final.haplotypes))]

            # k = comma seperated string of tids
            # v = the count
            for k, v in final.ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(final.target_idx_to_main_target[idx])

                #print(k, temp_main_targets)

                # loop through the targets and haplotypes to get the bits
                for main_target in temp_main_targets:
                    # main_target is not an index, but a value like 'ENMUST..001'

                    for i, hap in enumerate(final.haplotypes):

                        read_transcript = '{}_{}'.format(main_target, hap) # now 'ENMUST..001_A'

                        if len(final.haplotypes) == 1 and final.haplotypes[0] == '0':
                            read_transcript = main_target

                        # get the numerical tid corresponding to read_transcript
                        read_transcript_idx = str(alignment_file.gettid(read_transcript))

                        if read_transcript_idx in arr_target_idx:
                            #LOG.debug("{}\t{}\t{}".format(final.ec_idx[k], final.main_targets[main_target], i))

                            # main_targets[main_target] = idx of main target
                            # i = the haplotype
                            # ec_idx[k] = index of ec

                            #apm.set_value(final.main_targets[main_target], i, final.ec_idx[k], 1)

                            ec_arr[i].append(final.ec_idx[k])
                            target_arr[i].append(final.main_targets[main_target])
                            ec_arr_max = max(ec_arr_max, final.ec_idx[k])
                            target_arr_max = max(target_arr_max, final.main_targets[main_target])

            LOG.debug('APM shape={}'.format(new_shape))

            apm = APM(shape=new_shape,
                      haplotype_names=final.haplotypes,
                      locus_names=final.main_targets.keys(),
                      read_names=ec_ids,
                      sample_names=final.CRS.keys())

            for h in xrange(0, len(final.haplotypes)):
                d = np.ones(len(ec_arr[h]))
                apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])), shape=(len(final.ec), len(final.main_targets)))

            # now ER by UR
            # OLD: apm.count = final.ec.values()

            LOG.debug('Constructing CRS...')
            LOG.debug('CRS dimensions: {:,} x {:,}'.format(len(final.ec), len(final.CRS)))
            LOG.debug('Constructing CRS...')

            #print final.CRS_idx
            #print len(final.CRS)
            #print len(final.CRS_idx)

            npa = np.zeros((len(final.ec), len(final.CRS)))
            i = 0
            for eckey, crs in final.ec.iteritems():
                # k = commas seperated list (eckey)
                # v = dict of CRS and counts
                for crskey, crscount in crs.iteritems():
                    #print i, crskey, final.CRS_idx[crskey]
                    npa[i, final.CRS_idx[crskey]] = crs[crskey]
                i += 1

            LOG.info("NPA SUM: {:,}".format(npa.sum()))

            apm.count = npa #.reshape((len(final.ec), len(final.CRS)))

            LOG.info("APM Created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                utils.format_time(start_time, time.time())))

            temp_time = time.time()
            LOG.info("Flushing to disk...")
            apm.finalize()
            apm.save(output_filename, title='bam2ec')
            LOG.info("{} created in {}, total time: {}".format(output_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

        except Exception as e:
            LOG.fatal("ERROR: {}".format(str(e)))
            raise e
    else:
        try:
            temp_time = time.time()
            LOG.info("Generating BIN file...")

            f = open(output_filename, "wb")

            # version
            f.write(pack('<i', 2))
            LOG.info("VERSION: 2")

            # targets
            LOG.info("NUMBER OF TARGETS: {:,}".format(len(final.main_targets)))
            f.write(pack('<i', len(final.main_targets)))
            for main_target, idx in final.main_targets.iteritems():
                #LOG.debug("{:,}\t{}\t# {:,}".format(len(main_target), main_target, idx))
                f.write(pack('<i', len(main_target)))
                f.write(pack('<{}s'.format(len(main_target)), main_target))

            # haplotypes
            LOG.info("NUMBER OF HAPLOTYPES: {:,}".format(len(final.haplotypes)))
            f.write(pack('<i', len(final.haplotypes)))
            for idx, hap in enumerate(final.haplotypes):
                #LOG.debug("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
                f.write(pack('<i', len(hap)))
                f.write(pack('<{}s'.format(len(hap)), hap))


            #
            LOG.info("FILTERED CRS: {:,}".format(len(final.CRS)))
            f.write(pack('<i', len(final.CRS)))
            for CR, idx in final.CRS.iteritems():
                #LOG.debug("{:,}\t{}\t# {:,}".format(len(CR), CR, idx))
                f.write(pack('<i', len(CR)))
                f.write(pack('<{}s'.format(len(CR)), CR))

            # equivalence classes
            LOG.info("NUMBER OF EQUIVALANCE CLASSES: {:,}".format(len(final.ec)))
            LOG.info("NUMBER OF EQUIVALANCE CLASSES: {:,}".format(len(ec_totals)))
            f.write(pack('<i', len(ec_totals)))
            for idx, k in enumerate(ec_totals.keys()):
                # ec[k] is the count
                #LOG.debug("{:,}\t# {}\t{:,}".format(final.ec[k], k, idx))
                f.write(pack('<i', ec_totals[k]))

            LOG.debug("Determining mappings...")

            # equivalence class mappings
            counter = 0
            for k, v in final.ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(final.target_idx_to_main_target[idx])

                counter += len(temp_main_targets)

            LOG.info("NUMBER OF EQUIVALANCE CLASS MAPPINGS: {:,}".format(counter))
            f.write(pack('<i', counter))

            for k, v in final.ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(final.target_idx_to_main_target[idx])

                # loop through the haplotypes and targets to get the bits
                for main_target in temp_main_targets:
                    # main_target is not an index, but a value like 'ENMUST..001'

                    bits = []

                    # TODO: fix this, maybe set a flag
                    # checking to see if there is more than 1 haplotype and
                    # not the generic/default one of '0'
                    if len(final.haplotypes) == 1 and final.haplotypes[0] == '0':
                        read_transcript_idx = str(alignment_file.gettid(main_target))

                        if read_transcript_idx in arr_target_idx:
                            bits.append(1)
                        else:
                            bits.append(0)

                        #LOG.debug("{}\t{}\t{}\t# {}\t{}".format(final.ec_idx[k], final.main_targets[main_target], utils.list_to_int(bits), main_target, bits))
                        f.write(pack('<i', final.ec_idx[k]))
                        f.write(pack('<i', final.main_targets[main_target]))
                        f.write(pack('<i', utils.list_to_int(bits)))

                    else:
                        for hap in final.haplotypes:
                            read_transcript = '{}_{}'.format(main_target, hap) # now 'ENMUST..001_A'
                            read_transcript_idx = str(alignment_file.gettid(read_transcript))

                            if read_transcript_idx in arr_target_idx:
                                bits.append(1)
                            else:
                                bits.append(0)

                        #LOG.debug("{}\t{}\t{}\t# {}\t{}".format(final.ec_idx[k], final.main_targets[main_target], utils.list_to_int(bits), main_target, bits))
                        f.write(pack('<i', final.ec_idx[k]))
                        f.write(pack('<i', final.main_targets[main_target]))
                        f.write(pack('<i', utils.list_to_int(bits)))


            total = len(final.ec) * len(final.CRS)
            LOG.info("NUMBER OF RECORDS: {:,}".format(total))
            f.write(pack('<q', total))

            i = 0
            #for k, v in final.ec.iteritems():
            #    dump = [0] * len(final.CRS.keys())
            #    for idx, CRS in final.CRS.iteritems():
            #        if CRS in v:
            #            dump[idx] = v[CRS]

            for idx, CRS in final.CRS.iteritems():
                dump = [0] * len(final.ec.keys())
                for k, v in final.ec.iteritems():
                    if CRS in v:
                        dump[idx] = v[CRS]

                f.write(pack('<{}i'.format(len(dump)), *dump))
                #print(i)
                i += 1

            f.close()

            LOG.info("{} created in {}, total time: {}".format(output_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

        except Exception as e:
            LOG.error("Error: {}".format(str(e)))


def split_bam(filename, n, directory=None):
    """
    Split a BAM file into ``n`` files.

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
        aln = aln_file.next()
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
        f = open(filename, 'r')
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
            aln = aln_file.next()
            aln_qname = aln.query_name

            # added because some query_name had spaces
            i = aln_qname.find(' ')
            if i > 0:
                aln_qname = aln_qname[:i]

            aln_first_qname = aln_qname

            while aln_first_qname == aln_qname:

                virtual_offset = aln_file.tell()
                aln = aln_file.next()
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
        #print (params.input_file)

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
            for k, v in result.tid_ranges.iteritems():
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
