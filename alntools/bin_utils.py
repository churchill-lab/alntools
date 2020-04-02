# -*- coding: utf-8 -*-
from collections import OrderedDict
from six import iteritems
from struct import pack, unpack

import os
import time

import numpy as np

from scipy.sparse import coo_matrix, csr_matrix, csc_matrix

from .matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from . import bin_file
from . import utils

LOG = utils.get_logger()

try:
    xrange
except NameError:
    xrange = range


def get_ec_dict(self):
    ecs = OrderedDict()
    for idx, row in xrange(len(self.a_matrix.indptr)):
        ec_key = ','.join(map(str, self.a_matrix.getrow(idx).nonzero()[1]))
        ecs[ec_key] = len(ecs)


def load(ec_filename):
    with open(ec_filename, 'rb') as f:

        bin_format = unpack('<i', f.read(4))[0]

        if bin_format == 2:
        
            num_haps = unpack('<i', f.read(4))[0]
            hname = list()
            for hidx in range(num_haps):
                hname_len = unpack('<i', f.read(4))[0]
                hname.append(unpack('<{}s'.format(hname_len), f.read(hname_len))[0].decode('utf-8'))
            hname = np.array(hname)
        
            num_transcripts = unpack('<i', f.read(4))[0]
            transcript_lengths = np.zeros((num_transcripts, num_haps), dtype=float)
            tname = list()
            for tidx in range(num_transcripts):
                tname_len = unpack('<i', f.read(4))[0]
                tname.append(unpack('<{}s'.format(tname_len), f.read(tname_len))[0].decode('utf-8'))
                for hidx in range(num_haps):
                    transcript_lengths[tidx, hidx] = unpack('<i', f.read(4))[0]
            tname = np.array(tname)
        
            sname = list()
            num_samples = unpack('<i', f.read(4))[0]
            for sidx in range(num_samples):
                sname_len = unpack('<i', f.read(4))[0]
                sname.append(unpack('<{}s'.format(sname_len), f.read(sname_len))[0].decode('utf-8'))
            sname = np.array(sname)

            indptr_len = unpack('<i', f.read(4))[0]
            num_ecs = indptr_len - 1
            nnz = unpack('<i', f.read(4))[0]
            indptr_A = np.array(unpack('<{}i'.format(indptr_len), f.read(4*indptr_len)))
            indices_A = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
            data_A = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))

            indptr_len = unpack('<i', f.read(4))[0]
            nnz = unpack('<i', f.read(4))[0]
            indptr_N = np.array(unpack('<{}i'.format(indptr_len), f.read(4*indptr_len)))
            indices_N = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
            data_N = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))

            alnmat = APM()
            alnmat.shape = (num_transcripts, num_haps, num_ecs)
            alnmat.hname = hname
            alnmat.lname = tname
            alnmat.sname = sname
            alnmat.num_haplotypes = len(hname)
            alnmat.num_loci = len(tname)
            alnmat.num_reads = num_ecs
            alnmat.num_samples = len(sname)
            alnmat.lid = dict(zip(tname, np.arange(num_transcripts)))
            alnmat.sid = dict(zip(sname, np.arange(num_samples)))
            alnmat.lengths = transcript_lengths
            alnmat.finalized = False
            for hidx in range(num_haps-1):
                data_A, data_A_rem = np.divmod(data_A, 2)
                alnmat.data.append(csr_matrix((data_A_rem, indices_A, indptr_A), shape=(num_ecs, num_transcripts)))
            alnmat.data.append(csr_matrix((data_A, indices_A, indptr_A), shape=(num_ecs, num_transcripts)))
            for hidx in range(num_haps):
                alnmat.data[hidx].eliminate_zeros()
            alnmat.count = csc_matrix((data_N, indices_N, indptr_N), shape=(num_ecs, num_samples))
            if num_samples == 1:
                alnmat.count = alnmat.count.todense().A.flatten()
            alnmat.finalize()
            return alnmat

        elif bin_format == 1:
            raise NotImplementedError

        elif bin_format == 0:
            raise TypeError('Format 0 is not supported anymore.')


def save(ec_filename, sample, haplotypes, main_targets, main_target_lengths, alnmat, cntmat):
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
            f.write(pack('<{}s'.format(len(hap)), hap.encode('utf-8')))

        #
        # SECTION: TARGETS
        #     [# of TARGETS = T]
        #     [length TARGET 1 text][TARGET 1 text][HAP 1 length] ... [HAP H length]
        #     ...
        #     [length TARGET T text][TARGET T text][HAP 1 length] ... [HAP H length]
        #
        # Example:
        #     80000
        #     18 ENSMUST00000156068 234 235
        #     18 ENSMUST00000209341 1054 1054
        #     ...
        #     18 ENSMUST00000778019 1900 1890
        #

        LOG.info("NUMBER OF TARGETS: {:,}".format(len(main_targets)))
        f.write(pack('<i', len(main_targets)))
        for idx, main_target in enumerate(main_targets):
            f.write(pack('<i', len(main_target)))
            f.write(pack('<{}s'.format(len(main_target)), main_target.encode('utf-8')))
            for idx_hap, hap in enumerate(haplotypes):
                f.write(pack('<i', main_target_lengths[idx, idx_hap]))  # TODO: store as double (Fix emase-zero too)

        #
        # SECTION: CRS
        #     [# of CRS = 1]
        #     [length of CR 1 text][CR 1 text]
        #
        # Example:
        #     1
        #     8 SAMPLEID
        #

        LOG.info("NUMBER OF SAMPLES: 1")
        f.write(pack('<i', 1))
        f.write(pack('<i', len(sample)))
        f.write(pack('<{}s'.format(len(sample)), sample.encode('utf-8')))


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

        LOG.info("A MATRIX: INDPTR LENGTH {:,}".format(len(alnmat.indptr)))
        f.write(pack('<i', len(alnmat.indptr)))

        # NON ZEROS
        LOG.info("A MATRIX: NUMBER OF NON ZERO: {:,}".format(alnmat.nnz))
        f.write(pack('<i', alnmat.nnz))

        # ROW OFFSETS
        LOG.info("A MATRIX: LENGTH INDPTR: {:,}".format(len(alnmat.indptr)))
        f.write(pack('<{}i'.format(len(alnmat.indptr)), *alnmat.indptr))
        # LOG.error(alnmat.indptr)

        # COLUMNS
        LOG.info("A MATRIX: LENGTH INDICES: {:,}".format(len(alnmat.indices)))
        f.write(pack('<{}i'.format(len(alnmat.indices)), *alnmat.indices))
        # LOG.error(alnmat.indices)

        # DATA
        LOG.info("A MATRIX: LENGTH DATA: {:,}".format(len(alnmat.data)))
        f.write(pack('<{}i'.format(len(alnmat.data)), *alnmat.data))
        # LOG.error(alnmat.data)

        #
        # SECTION: "N" Matrix
        #
        # "N" Matrix format is EC (rows) by CRS (columns) with
        # each value being the EC count.
        #
        # Instead of storing a "dense" matrix, we store a "sparse"
        # matrix utilizing Compressed Sparse Column (CSC) format.
        #

        LOG.info("N MATRIX: NUMBER OF EQUIVALENCE CLASSES: {:,}".format(cntmat.shape[0]))
        LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(cntmat.indptr)))
        f.write(pack('<i', len(cntmat.indptr)))

        # NON ZEROS
        LOG.info("N MATRIX: NUMBER OF NON ZERO: {:,}".format(cntmat.nnz))
        f.write(pack('<i', cntmat.nnz))

        # ROW OFFSETS
        LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(cntmat.indptr)))
        f.write(pack('<{}i'.format(len(cntmat.indptr)), *cntmat.indptr))
        # LOG.error(cntmat.indptr)

        # COLUMNS
        LOG.info("N MATRIX: LENGTH INDICES: {:,}".format(len(cntmat.indices)))
        f.write(pack('<{}i'.format(len(cntmat.indices)), *cntmat.indices))
        # LOG.error(cntmat.indices)

        # DATA
        LOG.info("N MATRIX: LENGTH DATA: {:,}".format(len(cntmat.data)))
        f.write(pack('<{}i'.format(len(cntmat.data)), *cntmat.data))
        # LOG.error(cntmat.data)


def combine(ec_files, ec_out):
    start_time = time.time()

    # key = Haplotype name, value = Order inserted
    haplotypes = OrderedDict()
    haplotypes_idx = []

    # key = Transcript w/o Haplotype, value = Order inserted
    targets = OrderedDict()
    targets_idx = []

    targets_lengths = OrderedDict()

    # key = Sample, value = Order inserted
    samples = OrderedDict()
    samples_idx = []

    # key = ec, value = Order inserted
    ec_idx = {}

    ec_totals = OrderedDict()
    sample_totals = OrderedDict()

    ECS = OrderedDict()

    for ec_file_idx, ec_file_name in enumerate(ec_files):
        current_ECS = OrderedDict()

        #
        # parse file
        #
        ECF = bin_file.ECFile(ec_file_name)

        # since we are looping by rows, make the n matrix more efficient
        ECF.n_matrix = ECF.n_matrix.tocsr()

        ECF_num_haps = len(ECF.haplotypes_idx)

        temp_time = time.time()

        LOG.info('Looping through {:,} EC keys to generate dictionary'.format(ECF.a_matrix.shape[0]))

        for idx in xrange(ECF.a_matrix.shape[0]):

            if idx % 1000 == 0:
                LOG.info("idx: {}, time: {}, total time: {}".format(idx,
                                                                    utils.format_time(
                                                                        temp_time,
                                                                        time.time()),
                                                                    utils.format_time(start_time, time.time())))


            # get the row values
            a_row = ECF.a_matrix.getrow(idx)

            # only need the columns (targets) that have values
            a_col = list(a_row.nonzero()[1])

            #ec_key = ','.join(['{}:{}'.format(v, a_row[0, v]) for v in a_col])
            ec_key = ','.join(['{}:{}'.format(v, a_row[0, v]) for v in a_col])
            print('ec_key=', ec_key)

            # get the n matrix row (same ec)
            n_row = ECF.n_matrix.getrow(idx)

            # only need the columns (samples) that have values
            #print('n_row.nonzero()', n_row.nonzero())
            n_col = list(n_row.nonzero()[1])

            #print('idx=', idx)
            #print('n_row=', n_row)
            #print('n_col=', n_col)

            # EC[ec_key] = {SAMPLE: count...}
            #for v in n_col:
            #    print('n_row[0, v]=', n_row[0, v])
            #    print('ECF.samples_idx[v]=', ECF.samples_idx[v])

            current_ECS[ec_key] = {ECF.samples_idx[v]: n_row[0, v] for v in n_col}
            print(current_ECS[ec_key])

        #print('current_ECS=', current_ECS)

        LOG.info("Dictionary created {}, total time: {}".format(
                utils.format_time(temp_time, time.time()),
                utils.format_time(start_time, time.time())))

        temp_time = time.time()

        if ec_file_idx == 0:
            LOG.info("First file so generating target lengths, ec_totals, etc")

            ECS = current_ECS

            for (h, idx) in iteritems(ECF.haplotypes):
                haplotypes[h] = idx
                haplotypes_idx.append(h)

            for (s, idx) in iteritems(ECF.samples):
                samples[s] = idx
                samples_idx.append(s)

            for (t, idx) in iteritems(ECF.targets):
                targets[t] = idx
                targets_idx.append(t)

            for (t, h) in iteritems(ECF.targets_lengths):
                targets_lengths[t] = h

            #
            # we could probably do this more efficiently, but ...
            #
            for (eckey, crdict) in iteritems(ECS):
                ec_idx[eckey] = len(ec_idx)
                for (crkey, count) in iteritems(crdict):
                    if crkey in sample_totals:
                        sample_totals[crkey] += count
                    else:
                        sample_totals[crkey] = count

                    if eckey in ec_totals:
                        ec_totals[eckey] += count
                    else:
                        ec_totals[eckey] = count

            LOG.info("Done first file in {}, total time: {}".format(
                    utils.format_time(temp_time, time.time()),
                    utils.format_time(start_time, time.time())))

        else:
            #
            # we are assuming haplotypes and transcripts are the same
            #
            # haplotypes should be in same order
            #

            LOG.info("File number: {}, combining...".format(ec_file_idx))

            for (h, idx) in iteritems(ECF.haplotypes):
                if h not in haplotypes:
                    haplotypes[h] = len(haplotypes)
                    haplotypes_idx.append(h)

            for (s, idx) in iteritems(ECF.samples):
                if s not in samples:
                    samples[s] = len(samples)
                    samples_idx.append(s)

            for (t, idx) in iteritems(ECF.targets):
                if t not in targets:
                    targets[t] = len(targets)
                    targets_idx.append(t)

            for (t, h) in iteritems(ECF.targets_lengths):
                if t not in targets_lengths:
                    targets_lengths[t] = h

            for (eckey, crdict) in iteritems(current_ECS):
                for (crkey, count) in iteritems(crdict):
                    if crkey not in samples:
                        samples[crkey] = len(samples)

                    if crkey in sample_totals:
                        sample_totals[crkey] += count
                    else:
                        sample_totals[crkey] = count

                    if eckey in ec_totals:
                        ec_totals[eckey] += count
                    else:
                        ec_totals[eckey] = count

                    if eckey in ECS:
                        if crkey in ECS[eckey]:
                            ECS[eckey][crkey] += count
                        else:
                            ECS[eckey][crkey] = count
                    else:
                        ec_idx[eckey] = len(ec_idx)
                        ECS[eckey] = {crkey: count}

            LOG.info("Done file in {}, total time: {}".format(
                            utils.format_time(temp_time, time.time()),
                            utils.format_time(start_time, time.time())))

    #
    # create the binary file
    #

    #print('-------------------------------')
    #print('haplotypes = ', haplotypes)
    #print('targets = ', targets)
    #print('samples = ', samples)

    #print('ec_idx = ', ec_idx)
    #print('sample_totals = ', sample_totals)
    #print('ec_totals = ', ec_totals)

    try:
        temp_time = time.time()
        LOG.info('Constructing APM structure...')

        new_shape = (len(targets),
                     len(haplotypes),
                     len(ECS))

        LOG.debug('Shape={}'.format(new_shape))

        # final.ec.values -> the number of times this equivalence class has appeared

        ec_ids = [x for x in xrange(0, len(ECS))]
        ec_arr = [[] for _ in xrange(0, len(haplotypes))]
        target_arr = [[] for _ in xrange(0, len(haplotypes))]

        #print('ec_ids=', ec_ids)
        #print('ec_arr=', ec_arr)
        #print('target_arr=', target_arr)


        indptr = [0]
        indices = []
        data = []

        #print('ec_idx=', ec_idx)
        #print('targets=', targets)

        # k = comma seperated string of tids:count
        # v = dict of samples and counts
        for (eckey, crdict) in iteritems(ECS):
            #print('')
            #print('eckey=', eckey)
            #print('crdict=', crdict)
            target_info = eckey.split(',')
            #print('target_info=', target_info)

            for idx, target in enumerate(target_info):
                elems = target.split(':')

                target = elems[0]
                haplotype = elems[1]
                #print('target, haplotype=', target, haplotype)

                haps = utils.int_to_list(int(haplotype), ECF_num_haps)
                #print('haps=', haps)

                for h_idx, h in enumerate(haps):
                    if h != 0:
                        #print('h_idx, h = ', h_idx, h)

                        ec_arr[h_idx].append(ec_idx[eckey])
                        target_arr[h_idx].append(int(target))

            # construct "N" matrix elements
            ti = sorted(crdict.keys(), key=lambda i: samples[i])

            a = 0
            for crskey in ti:
                #print('crs_key=', crskey)
                col = samples[crskey]
                #print('col=', col)
                indices.append(col)
                #print('crdict[crskey]=', crdict[crskey])
                data.append(crdict[crskey])
                a += 1

            indptr.append(indptr[-1] + a)

        apm = APM(shape=new_shape,
                  haplotype_names=haplotypes,
                  locus_names=targets.keys(),
                  read_names=ec_ids,
                  sample_names=samples.keys())

        for h in xrange(0, len(haplotypes)):
            d = np.ones(len(ec_arr[h]), dtype=np.int32)
            apm.data[h] = coo_matrix((d, (ec_arr[h], target_arr[h])),
                                     shape=(len(ECS), len(targets)))

        LOG.debug('Constructing CRS...')
        LOG.debug('CRS dimensions: {:,} x {:,}'.format(len(ECS), len(samples)))

        LOG.info('len data={}'.format(len(data)))
        LOG.info('data={}'.format(data[-10:]))
        LOG.info('len indices={}'.format(len(indices)))
        LOG.info('indices={}'.format(indices[-10:]))
        LOG.info('len indptr={}'.format(len(indptr)))
        LOG.info('indptr={}'.format(indptr[-10:]))

        npa = csr_matrix((np.array(data, dtype=np.int32),
                          np.array(indices, dtype=np.int32),
                          np.array(indptr, dtype=np.int32)),
                          shape=(len(ECS), len(samples)))

        LOG.info("NPA SUM: {:,}".format(npa.sum()))

        apm.count = npa.tocsc()

        LOG.info("APM Created in {}, total time: {}".format(
            utils.format_time(temp_time, time.time()),
            utils.format_time(start_time, time.time())))

        '''

        if emase_filename:
            LOG.info("Flushing to disk...")

            try:
                os.remove(emase_filename)
            except OSError:
                pass

            temp_time = time.time()
            apm.finalize()
            apm.save(emase_filename, title='Multisample APM',
                     incidence_only=False)
            LOG.info("{} created in {}, total time: {}".format(emase_filename,
                                                               utils.format_time(
                                                                   temp_time,
                                                                   time.time()),
                                                               utils.format_time(
                                                                   start_time,
                                                                   time.time())))
        '''

        if ec_out:
            LOG.debug("Creating summary matrix...")

            try:
                os.remove(ec_out)
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

            LOG.info("Matrix created in {}, total time: {}".format(
                utils.format_time(temp_time, time.time()),
                utils.format_time(start_time, time.time())))

            temp_time = time.time()
            LOG.info("Generating BIN file...")

            with open(ec_out, 'wb') as f:
                # FORMAT
                f.write(pack('<i', 2))
                LOG.info("FORMAT: 2")

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

                LOG.info("NUMBER OF HAPLOTYPES: {:,}".format(len(haplotypes)))
                f.write(pack('<i', len(haplotypes)))
                for (hap, idx) in iteritems(haplotypes):
                    # LOG.debug("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
                    f.write(pack('<i', len(hap)))
                    f.write(pack('<{}s'.format(len(hap)), hap))

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

                LOG.info("NUMBER OF TARGETS: {:,}".format(len(targets)))
                f.write(pack('<i', len(targets)))
                for (main_target, idx) in iteritems(targets):
                    f.write(pack('<i', len(main_target)))
                    f.write(pack('<{}s'.format(len(main_target)), main_target))

                    # lengths = []

                    for (hap, idx_hap) in iteritems(haplotypes):
                        length = targets_lengths[main_target][hap]
                        f.write(pack('<i', length))
                        # lengths.append(str(length))

                    # LOG.debug("#{:,} --> {:,}\t{}\t{}\t".format(idx, len(main_target), main_target, '\t'.join(lengths)))

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

                LOG.info("FILTERED CRS: {:,}".format(len(samples)))
                f.write(pack('<i', len(samples)))
                for (sample, idx) in iteritems(samples):
                    # LOG.debug("{:,}\t{}\t# {:,}".format(len(CR), CR, idx))
                    f.write(pack('<i', len(sample)))
                    f.write(pack('<{}s'.format(len(sample)), sample))

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

                LOG.info("Determining mappings...")

                num_mappings = summat.nnz
                summat = summat.tocsr()

                LOG.info("A MATRIX: INDPTR LENGTH {:,}".format(len(summat.indptr)))
                f.write(pack('<i', len(summat.indptr)))

                # NON ZEROS
                LOG.info("A MATRIX: NUMBER OF NON ZERO: {:,}".format(num_mappings))
                f.write(pack('<i', num_mappings))

                # ROW OFFSETS
                LOG.info("A MATRIX: LENGTH INDPTR: {:,}".format(len(summat.indptr)))
                f.write(pack('<{}i'.format(len(summat.indptr)), *summat.indptr))
                LOG.debug(summat.indptr)

                # COLUMNS
                LOG.info("A MATRIX: LENGTH INDICES: {:,}".format(len(summat.indices)))
                f.write(pack('<{}i'.format(len(summat.indices)), *summat.indices))
                LOG.debug(summat.indices)

                # DATA
                LOG.info("A MATRIX: LENGTH DATA: {:,}".format(len(summat.data)))
                f.write(pack('<{}i'.format(len(summat.data)), *summat.data))
                LOG.debug(summat.data)

                #
                # SECTION: "N" Matrix
                #
                # "N" Matrix format is EC (rows) by CRS (columns) with
                # each value being the EC count.
                #
                # Instead of storing a "dense" matrix, we store a "sparse"
                # matrix utilizing Compressed Sparse Column (CSC) format.
                #

                LOG.info("N MATRIX: NUMBER OF EQUIVALENCE CLASSES: {:,}".format(len(ECS)))
                LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<i', len(apm.count.indptr)))

                # NON ZEROS
                LOG.info("N MATRIX: NUMBER OF NON ZERO: {:,}".format(apm.count.nnz))
                f.write(pack('<i', apm.count.nnz))

                # ROW OFFSETS
                LOG.info("N MATRIX: LENGTH INDPTR: {:,}".format(len(apm.count.indptr)))
                f.write(pack('<{}i'.format(len(apm.count.indptr)), *apm.count.indptr))
                LOG.debug(apm.count.indptr)

                # COLUMNS
                LOG.info("N MATRIX: LENGTH INDICES: {:,}".format(len(apm.count.indices)))
                f.write(pack('<{}i'.format(len(apm.count.indices)), *apm.count.indices))
                LOG.debug(apm.count.indices)

                # DATA
                LOG.info("N MATRIX: LENGTH DATA: {:,}".format(len(apm.count.data)))
                f.write(pack('<{}i'.format(len(apm.count.data)), *apm.count.data))
                LOG.debug(apm.count.data)

            LOG.info("{} created in {}, total time: {}".format(ec_out,
                                                               utils.format_time(temp_time, time.time()),
                                                               utils.format_time(start_time,time.time())))
    except KeyboardInterrupt as e:
        LOG.error("Error: {}".format(str(e)))


def bin2apm(bin_filename, apm_filename):
    try:
        start_time = time.time()

        ECF = bin_file.ECFile(bin_filename)
        apm = ECF.toAPM()

        LOG.info("Flushing to disk...")

        try:
            os.remove(apm_filename)
        except OSError:
            pass

        temp_time = time.time()
        apm.finalize()
        apm.save(apm_filename, title='Multisample APM', incidence_only=False)
        LOG.info("{} created in {}, total time: {}".format(apm_filename,
                                                           utils.format_time(temp_time, time.time()),
                                                           utils.format_time(start_time, time.time())))

    except Exception as e:
        LOG.error("Error: {}".format(str(e)))
