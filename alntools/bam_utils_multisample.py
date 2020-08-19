# -*- coding: utf-8 -*-
from collections import OrderedDict, namedtuple, defaultdict
from struct import pack

import glob
import multiprocessing
import os
import sys
import time

from six import iteritems

from scipy.sparse import coo_matrix, lil_matrix, csr_matrix

import numpy as np
import pysam

from . import utils
from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

try:
    xrange
except NameError:
    xrange = range


LOG = utils.get_logger()
BAM_HEADER = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00'
BAM_EOF = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'

parse_fields = ["header_size", 
                "begin_read_offset", 
                "begin_read_size", 
                "file_offset", 
                "file_bytes", 
                "end_read_offset",
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
    slots = ['input_file', 'temp_dir', 'process_id', 'data', 'track_ranges']

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
    slots = ['valid_alignments', 'all_alignments', 'ec', 'init', 'tid_ranges', 'CRS']

    def __init__(self):
        self.valid_alignments = None
        self.all_alignments = None
        self.ec = None
        self.init = False
        self.tid_ranges = None
        self.CRS = None


class RangeParams(object):
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
    for (k, v) in iteritems(d):
        if isinstance(v, dict):
            d[k] = ddict2dict(v)
    return dict(d)


def process_convert_bam(cp):
    LOG.debug('Process {}: Input File: {}'.format(cp.process_id, cp.input_file))

    validate_bam(cp.input_file)

    main_targets = OrderedDict()

    # reference_id: the reference sequence number as defined in the header
    # reference_name: name (None if no AlignmentFile is associated)

    # ec = equivalence class
    #      the KEY is a comma separated string of reference_ids
    #      the VALUE is the number of times this equivalence class has appeared
    ec = defaultdict(lambda: defaultdict(int))

    # unique reads
    read_IDs = {}

    # times encountering new read id
    read_id_switch_counter = 0
    same_read_target_counter = 0

    ranges = {}

    all_alignments = 0
    valid_alignments = 0
    ec_key = None
    temp_name = os.path.join(cp.temp_dir, '_bam2ec.')
    generic_haplotype = ''

    reference_id = None
    reference_ids = []

    try:
        alignment_file = pysam.AlignmentFile(cp.input_file)

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

            if cp.track_ranges:
                min_max = ranges.get(reference_id, (100000000000, -1))
                n = min(min_max[0], alignment.reference_start)
                x = max(min_max[1], alignment.reference_start)
                ranges[reference_id] = (n,x)

            # query_name = Column 1 from file, the Query template NAME

            if query_name is None:
                query_name = alignment.query_name

                i = query_name.find(' ')
                if i > 0:
                    query_name = query_name[:i]

            try:
                read_IDs[query_name] += 1
            except KeyError:
                read_IDs[query_name] = 1

            # query_name = Column 1 from file, the Query template NAME
            bam_tags = query_name.split('|||')
            #orig_query_name = bam_tags[0]
            #bam_tag_CR = bam_tags[2]
            #bam_tag_CY = bam_tags[4]
            #bam_tag_UR = bam_tags[6]
            #bam_tag_UY = bam_tags[8]
            #bam_tag_BC = bam_tags[10]
            #bam_tag_QT = bam_tags[12]

            # NEW CID
            bam_tag_CID = bam_tags[14]

            # set tags and dump new alignment
            alignment_query_name = alignment.query_name
            i = alignment_query_name.find(' ')
            if i > 0:
                alignment_query_name = alignment_query_name[:i]

            if query_name != alignment_query_name:
                ec_key = ','.join(sorted(reference_ids))
                ec[ec_key][bam_tag_CID] += 1

                query_name = alignment.query_name

                reference_ids = [reference_id]
                read_id_switch_counter += 1
            else:
                if reference_id not in reference_ids:
                    reference_ids.append(reference_id)
                else:
                    same_read_target_counter += 1

            if all_alignments % 100000 == 0:
                LOG.debug("Process {}: Done processing file {} ({:,} valid alignments processed out of {:,}, with {:,} equivalence classes)".format(
                          cp.process_id, cp.input_file, valid_alignments, all_alignments, len(ec)))

    except StopIteration:
        LOG.info("Process {}: Done processing file {} ({:,} valid alignments out of {:,}, with {:,} equivalence classes)".format(
                 cp.process_id, cp.input_file, valid_alignments, all_alignments, len(ec)))

    LOG.debug("Process {}: Number of Unique Reads: {:,}".format(cp.process_id, len(read_IDs)))
    LOG.debug("Process {}: Number of Reads/Target Duplications: {:,}".format(cp.process_id, same_read_target_counter))
    LOG.debug("Process {}: Number of Equivalence Classes: {:,}".format(cp.process_id, len(ec)))

    ret = ConvertResults()
    ret.valid_alignments = valid_alignments
    ret.all_alignments = all_alignments
    ret.ec = ddict2dict(ec)
    ret.read_IDs = read_IDs
    ret.tid_ranges = ranges

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

def convert(bam_filename, ec_filename, emase_filename, num_chunks=0, minimum_count=1, number_processes=-1, temp_dir=None, range_filename=None, target_filename=None):
    """
    """
    LOG.debug('Parameters')
    LOG.debug('-------------------------------------------')
    LOG.debug('BAM file: {}'.format(bam_filename))
    LOG.debug('EC file: {}'.format(ec_filename))
    LOG.debug('EMASE file: {}'.format(emase_filename))
    LOG.debug('chunks: {}'.format(num_chunks))
    LOG.debug('minimum count: {}'.format(minimum_count))
    LOG.debug('processes: {}'.format(number_processes))
    LOG.debug('Temp directory: {}'.format(temp_dir))
    LOG.debug('Range file: {}'.format(range_filename))
    LOG.debug('Target file: {}'.format(target_filename))
    LOG.debug('-------------------------------------------')

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
        if emase_filename:
            temp_dir = os.path.dirname(emase_filename)
        if ec_filename:
            temp_dir = os.path.dirname(ec_filename)

    LOG.info("Parsing the header of {}...".format(bam_files[0]))
    temp_time = time.time()

    alignment_file = pysam.AlignmentFile(bam_files[0])
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
    alignment_file = pysam.AlignmentFile(bam_files[0])

    #
    # Due to pysam limitation (at least on 0.10.0) we have to LOOP THROUGH TWICE
    # because referencing both alignment_file.lengths and
    # because referencing both alignment_file.references and
    #
    for idx, length in enumerate(alignment_file.lengths):
        reference_sequence_name = all_targets_list[idx]

        idx_underscore = reference_sequence_name.rfind('_')

        if idx_underscore > 0:
            target = reference_sequence_name[:idx_underscore]
            haplotype = reference_sequence_name[idx_underscore + 1:]
        else:
            target = reference_sequence_name
            haplotype = ''

        main_target_lengths[main_targets[target], haplotypes_idx[haplotype]] = length

    LOG.info("File parsed in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                        utils.format_time(start_time, time.time())))

    all_params = []

    pid = 0
    for bam_file in bam_files:
        params = ConvertParams()
        params.input_file = bam_file
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
    final.read_IDs = {}
    final.tid_ranges = {}

    ec_totals = OrderedDict()
    cr_totals = OrderedDict()

    LOG.info("Starting {} processes ...".format(num_processes))

    temp_time = time.time()
    args = zip(all_params)
    pool = multiprocessing.Pool(num_processes)

    ec_idx = {}
    CRS = OrderedDict()

    # parse results
    for idx, result in enumerate(pool.imap(wrapper_convert, args)):
        LOG.info("Process {} done out of {}. Combining results...".format(idx + 1, len(all_params)))
        temp_time = time.time()

        if not final.init:
            final = result
            final.init = True
            CRS = OrderedDict()
            new_ecs = OrderedDict()

            for (eckey, crdict) in iteritems(final.ec):
                new_ecs[eckey] = crdict
                ec_idx[eckey] = len(ec_idx)
                for (crkey, count) in iteritems(crdict):
                    if crkey not in CRS:
                        CRS[crkey] = len(CRS)

                    if crkey in cr_totals:
                        cr_totals[crkey] += count
                    else:
                        cr_totals[crkey] = count

                    if eckey in ec_totals:
                        ec_totals[eckey] += count
                    else:
                        ec_totals[eckey] = count

            final.ec = new_ecs

        else:
            # combine ec and URS
            LOG.debug("CHUNK {}: # Result Equivalence Classes: {:,}".format(idx, len(result.ec)))

            for (eckey, crdict) in iteritems(result.ec):
                for (crkey, count) in iteritems(crdict):
                    CRS[crkey] = 1

                    if crkey not in CRS:
                        CRS[crkey] = len(CRS)

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
                    else:
                        ec_idx[eckey] = len(ec_idx)
                        final.ec[eckey] = {crkey: count}

            LOG.debug("CHUNK {}: # Total Equivalence Classes: {:,}".format(idx, len(final.ec)))
            LOG.debug("CHUNK {}: # Total CRs: {:,}".format(idx, len(CRS)))

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

    LOG.info("Number of reads processed: {:,}".format(len(final.read_IDs)))
    # LOG.info("Number of total alignments: {:,}".format(final.all_alignments))
    LOG.info("Number of alignments: {:,}".format(final.valid_alignments))
    LOG.info("Number of main targets: {:,}".format(len(main_targets)))
    LOG.info("Number of haplotypes: {:,}".format(len(haplotypes)))
    LOG.info("Number of ECs: {:,}".format(len(final.ec)))
    LOG.info("Number of cells: {:,}".format(len(CRS)))

    LOG.debug("Number of ECs (ec_idx): {:,}".format(len(ec_idx)))
    LOG.debug("ECs Max Index: {:,}".format(max(ec_idx.values())))

    # keep all the cells (not recommended)
    if minimum_count <= 0:
        minimum_count = 1

    LOG.info("Minimum UMI count (or read count) per cell: {:,}".format(minimum_count))
    LOG.info("Filtering cells...")
    
    # find the new CRS
    CRS = OrderedDict()

    for (CR, CR_total) in iteritems(cr_totals):
        if CR_total >= minimum_count:
            if CR not in CRS:
                CRS[CR] = len(CRS)

    # remove invalid CRS from ECs
    new_ecs = OrderedDict()
    new_ec_idx = {}
    new_ec_totals = {}
    
    # loop through ecs
    for (eckey, crs) in iteritems(final.ec):
        # potential new ec

        ec = {}
        total = 0
        # loop through valid CRS and and if valid
        for (crkey, crcount) in iteritems(crs):

            if crkey in CRS:
                ec[crkey] = crs[crkey]
                total += crcount

        # only add to new ecs if there is anything
        if len(ec) > 0:
            new_ecs[eckey] = ec
            new_ec_idx[eckey] = len(new_ec_idx)
            new_ec_totals[eckey] = total

    #ec_totals = new_ec_totals
    final.ec = new_ecs
    ec_idx = new_ec_idx

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

    # LOG.info("# Valid Alignments: {:,}".format(final.valid_alignments))
    LOG.info("Number of ECs after filtering : {:,}".format(len(final.ec)))
    LOG.info("Number of cells after filtering: {:,}".format(len(CRS)))
    LOG.debug("Number of ECs (ec_idx) after filtering: {:,}".format(len(ec_idx)))
    LOG.debug("ECs Max Index after filtering: {:,}".format(max(ec_idx.values())))


    try:
        temp_time = time.time()
        # LOG.info('Constructing APM structure...')
        LOG.info('Constructing A matrix (alignment incidence matrix)...')

        new_shape = (len(main_targets),
                     len(haplotypes),
                     len(final.ec))

        LOG.debug('Shape={}'.format(new_shape))

        # final.ec.values -> the number of times this equivalence class has appeared

        ec_ids = [x for x in xrange(0, len(final.ec))]
        ec_arr = [[] for _ in xrange(0, len(haplotypes))]
        target_arr = [[] for _ in xrange(0, len(haplotypes))]

        indptr = [0]
        indices = []
        data = []

        # k = comma seperated string of tids
        # v = CRS and counts
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
                        #LOG.debug("{}\t{}\t{}".format(ec_idx[k], main_targets[main_target], i))

                        # main_targets[main_target] = idx of main target
                        # i = the haplotype
                        # ec_idx[k] = index of ec

                        #apm.set_value(main_targets[main_target], i, ec_idx[k], 1)

                        ec_arr[i].append(ec_idx[k])
                        target_arr[i].append(main_targets[main_target])

            # construct "N" matrix elements
            ti = sorted(v.keys(), key=lambda i: CRS[i])

            a = 0
            for crskey in ti:
                col = CRS[crskey]
                indices.append(col)
                data.append(v[crskey])
                a += 1

            indptr.append(indptr[-1] + a)

        apm = APM(shape=new_shape,
                  haplotype_names=haplotypes,
                  locus_names=main_targets.keys(),
                  read_names=ec_ids.astype(str),
                  sample_names=CRS.keys())

        for h in xrange(0, len(haplotypes)):
            d = np.ones(len(ec_arr[h]), dtype=np.int32)
            apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])), shape=(len(final.ec), len(main_targets)))

        LOG.info('Constructing N matrix (EC count matrix)...')
        LOG.debug('N matrix dimensions: {:,} x {:,}'.format(len(final.ec), len(CRS)))

        '''
        npa = lil_matrix((len(final.ec), len(CRS)), dtype=np.int32)
        i = 0
        for eckey, crs in final.ec.iteritems():
            # eckey = commas seperated list
            # crs = dict of CRS and counts
            for crskey, crscount in crs.iteritems():
                npa[i, CRS_idx[crskey]] = crs[crskey]
            i += 1
        '''
        LOG.debug('len data={}'.format(len(data)))
        LOG.debug('data={}'.format(data[-10:]))
        LOG.debug('len indices={}'.format(len(indices)))
        LOG.debug('indices={}'.format(indices[-10:]))
        LOG.debug('len indptr={}'.format(len(indptr)))
        LOG.debug('indptr={}'.format(indptr[-10:]))

        npa = csr_matrix((np.array(data, dtype=np.int32),
                          np.array(indices, dtype=np.int32),
                          np.array(indptr, dtype=np.int32)),
                         shape=(len(final.ec), len(CRS)))

        LOG.info("Number of UMIs in total: {:,}".format(npa.sum()))
        LOG.debug("NPA SUM: {:,}".format(npa.sum()))

        apm.count = npa.tocsc()

        apm.finalize()
        LOG.debug("APM Created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                             utils.format_time(start_time, time.time())))

        if emase_filename:
            LOG.info("Saving to {}...".format(emase_filename))

            try:
                os.remove(emase_filename)
            except OSError:
                pass

            temp_time = time.time()
            apm.save(emase_filename, title='Multisample APM', incidence_only=False)
            LOG.info("{} created in {}, total time: {}".format(emase_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

        #
        # TODO: Refactor using scsave method
        #
        if ec_filename:
            LOG.info("Saving to {}...".format(ec_filename))

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

            LOG.debug("Matrix created in {}, total time: {}".format(utils.format_time(temp_time, time.time()),
                                                                    utils.format_time(start_time, time.time())))

            temp_time = time.time()

            with open(ec_filename, 'wb') as f:
                # FORMAT
                f.write(pack('<i', 2))
                LOG.debug("EC FORMAT: 2")

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

                LOG.debug("Number of haplotypes: {:,}".format(len(haplotypes)))
                f.write(pack('<i', len(haplotypes)))
                for idx, hap in enumerate(haplotypes):
                    # LOG.debug("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
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

                LOG.debug("Number of reference targets: {:,}".format(len(main_targets)))
                f.write(pack('<i', len(main_targets)))
                for (main_target, idx) in iteritems(main_targets):
                    f.write(pack('<i', len(main_target)))
                    f.write(pack('<{}s'.format(len(main_target)), main_target))

                    #lengths = []

                    for idx_hap, hap in enumerate(haplotypes):
                        length = main_target_lengths[idx, idx_hap]
                        f.write(pack('<i', length))
                        #lengths.append(str(length))

                    #LOG.debug("#{:,} --> {:,}\t{}\t{}\t".format(idx, len(main_target), main_target, '\t'.join(lengths)))

                #
                # SECTION: CRS
                #     [# of CRS = C]
                #     [length of CR 1 text][CR 1 text]
                #     ...
                #     [length of CR C text][CR C text]
                #
                # Example:
                #     3
                #     16 TCGGTAAAGCCGTCGT
                #     16 GGAACTTAGCCGATTT
                #     16 TAGTGGTAGAGGTAGA
                #

                LOG.debug("Number of cells: {:,}".format(len(CRS)))
                f.write(pack('<i', len(CRS)))
                for (CR, idx) in iteritems(CRS):
                    #LOG.debug("{:,}\t{}\t# {:,}".format(len(CR), CR, idx))
                    f.write(pack('<i', len(CR)))
                    f.write(pack('<{}s'.format(len(CR)), CR))

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

                LOG.info("Saving A matrix...")

                num_mappings = summat.nnz
                summat = summat.tocsr()

                LOG.debug("A MATRIX: INDPTR LENGTH {:,}".format(len(summat.indptr)))
                f.write(pack('<i', len(summat.indptr)))

                # NON ZEROS
                LOG.debug("A MATRIX: NUMBER OF NON ZERO: {:,}".format(num_mappings))
                f.write(pack('<i', num_mappings))

                # ROW OFFSETS
                LOG.debug("A MATRIX: LENGTH INDPTR: {:,}".format(len(summat.indptr)))
                f.write(pack('<{}i'.format(len(summat.indptr)), *summat.indptr))
                # LOG.error(summat.indptr)

                # COLUMNS
                LOG.debug("A MATRIX: LENGTH INDICES: {:,}".format(len(summat.indices)))
                f.write(pack('<{}i'.format(len(summat.indices)), *summat.indices))
                # LOG.error(summat.indices)

                # DATA
                LOG.debug("A MATRIX: LENGTH DATA: {:,}".format(len(summat.data)))
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

                LOG.info("Saving N matrix...")

                LOG.debug("N MATRIX: NUMBER OF EQUIVALENCE CLASSES: {:,}".format(len(final.ec)))
                LOG.debug("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<i', len(apm.count.indptr)))

                # NON ZEROS
                LOG.debug("N MATRIX: NUMBER OF NON ZERO: {:,}".format(apm.count.nnz))
                f.write(pack('<i', apm.count.nnz))

                # ROW OFFSETS
                LOG.debug("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<{}i'.format(len(apm.count.indptr)), *apm.count.indptr))
                # LOG.error(apm.count.indptr)

                # COLUMNS
                LOG.debug("N MATRIX: LENGTH INDICES: {:,}".format(len(apm.count.indices)))
                f.write(pack('<{}i'.format(len(apm.count.indices)), *apm.count.indices))
                # LOG.error(apm.count.indices)

                # DATA
                LOG.debug("N MATRIX: LENGTH DATA: {:,}".format(len(apm.count.data)))
                f.write(pack('<{}i'.format(len(apm.count.data)), *apm.count.data))
                # LOG.error(apm.count.data)

            LOG.info("{} created in {}, total time: {}".format(ec_filename,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time, time.time())))

    except KeyboardInterrupt as e:
        LOG.error("Error: {}".format(str(e)))
