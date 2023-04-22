from collections import OrderedDict
import os
import time

from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
import numpy as np

from alntools.matrix.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM
from alntools import utils
from alntools import bin_utils

LOG = utils.get_logger()


def parse_salmon_ec(salmon_dir, target_filename=None):
    LOG.debug("-------------------------------------------")
    LOG.debug("Parameters:")
    LOG.debug(f"  SALMON directory: {salmon_dir}")
    LOG.debug(f"  Target file: {target_filename}")
    LOG.debug("-------------------------------------------")

    if target_filename is not None:
        LOG.info(
            f"Reading in the reference transcript IDs from {target_filename}"
        )
        transcripts_full = np.loadtxt(
            target_filename, dtype=str, delimiter="\t", usecols=(0,)
        )

    transcript_idx = OrderedDict()
    haplotype_idx = OrderedDict()
    tidx2coord = OrderedDict()
    targets = list()
    quant_file = os.path.join(salmon_dir, "quant.sf")
    salmon_file = os.path.join(salmon_dir, "aux_info", "eq_classes.txt")
    with open(salmon_file) as fh:
        LOG.info(f"Parsing {salmon_file}")
        num_targets = int(fh.readline())
        num_ec = int(fh.readline())
        tidx = 0
        hidx = 0
        for _ in range(num_targets):
            curline = fh.readline()
            tkey = curline.rstrip()
            targets.append(tkey)
            tgt, hap = tkey.split("_")
            if tgt not in transcript_idx:
                transcript_idx[tgt] = tidx
                tidx += 1
            if hap not in haplotype_idx:
                haplotype_idx[hap] = hidx
                hidx += 1
        if target_filename is not None:
            for tgt in transcripts_full:
                if not tgt in transcript_idx:
                    transcript_idx[tgt] = tidx
                    tidx += 1
        for tidx, tkey in enumerate(targets):
            tgt, hap = tkey.split("_")
            tidx2coord[tidx] = [transcript_idx[tgt], haplotype_idx[hap]]
        num_haps = len(haplotype_idx)
        num_transcripts = len(transcript_idx)

        LOG.info(
            f"Reading in the effective transcript lengths from {quant_file}"
        )
        targets = dict(zip(targets, np.arange(num_targets)))
        with open(quant_file) as qfh:
            qfh.readline()
            for curline in qfh:
                item = curline.rstrip().split("\t")
                tidx2coord[targets[item[0]]].append(float(item[2]))

        LOG.info("Creating EC alignment incidence matrix")
        rowid = list()
        colid = list()
        for h in range(num_haps):
            rowid.append(list())
            colid.append(list())

        ecidx = 0
        ec_counts = list()
        for curline in fh:
            item = curline.rstrip().split("\t")
            ec_counts.append(int(item[-1]))
            for ecitem in item[1:-1]:
                tidx, hidx, _ = tidx2coord[int(ecitem)]
                rowid[hidx].append(ecidx)
                colid[hidx].append(tidx)
            ecidx += 1

        data = list()
        for h in range(num_haps):
            data.append(
                coo_matrix(
                    (np.ones(len(rowid[h])), (rowid[h], colid[h])),
                    shape=(num_ec, num_transcripts),
                    dtype=int,
                )
            )
        for h in range(num_haps):
            data[h] = data[h].tocsr()

        alnmat = data[0]
        for h in range(1, num_haps):
            alnmat = alnmat + ((2**h) * data[h])

        LOG.info("Creating EC count matrix")
        ec_counts = np.array(ec_counts)
        cntmat = csr_matrix(
            ec_counts
        ).transpose()  # Note: This is eventually a csc_matrix

        LOG.info("Creating objects of meta data")
        transcript_lengths = np.zeros((num_transcripts, num_haps), dtype=int)
        for _, tinfo in tidx2coord.items():
            transcript_lengths[tinfo[0], tinfo[1]] = tinfo[2]
        # transcript_lengths[transcript_lengths < 1] = 1.0
        haplotypes = np.array(list(haplotype_idx.keys()))
        transcripts = np.array(list(transcript_idx.keys()))

    LOG.info("Parsing complete")
    return transcripts, haplotypes, alnmat, cntmat, transcript_lengths


# def parse_salmon_ec(salmon_dir, target_filename=None):
def parse_salmon_ec2(salmon_dir, target_filename=None):
    LOG.debug("-------------------------------------------")
    LOG.debug("Parameters:")
    LOG.debug(f"  SALMON directory: {salmon_dir}")
    LOG.debug(f"  Target file: {target_filename}")
    LOG.debug("-------------------------------------------")

    alnmat = APM()

    if target_filename is not None:
        LOG.info(
            f"Reading in the reference transcript IDs from {target_filename}"
        )
        transcripts_full = np.loadtxt(
            target_filename, dtype=str, delimiter="\t", usecols=(0,)
        )

    transcript_idx = OrderedDict()
    haplotype_idx = OrderedDict()
    tidx2coord = OrderedDict()
    targets = list()
    quant_file = os.path.join(salmon_dir, "quant.sf")
    salmon_file = os.path.join(salmon_dir, "aux_info", "eq_classes.txt")
    with open(salmon_file) as fh:
        LOG.info(f"Parsing {salmon_file}")
        alnmat.num_loci = int(fh.readline())
        alnmat.num_reads = int(fh.readline())
        tidx = 0
        hidx = 0
        for _ in range(alnmat.num_loci):
            curline = fh.readline()
            tkey = curline.rstrip()
            targets.append(tkey)
            tgt, hap = tkey.split("_")
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
            tgt, hap = tkey.split("_")
            tidx2coord[tidx] = [transcript_idx[tgt], haplotype_idx[hap]]

        alnmat.num_haplotypes = len(haplotype_idx)
        alnmat.num_samples = 1

        LOG.info(
            f"Reading in the effective transcript lengths from {quant_file}"
        )
        targets = dict(zip(targets, np.arange(alnmat.num_loci)))
        with open(quant_file) as qfh:
            qfh.readline()
            for curline in qfh:
                item = curline.rstrip().split("\t")
                tidx2coord[targets[item[0]]].append(float(item[2]))

        LOG.info("Creating EC alignment incidence matrix")
        rowid = list()
        colid = list()
        for h in range(alnmat.num_haplotypes):
            rowid.append(list())
            colid.append(list())

        ecidx = 0
        ec_counts = list()
        for curline in fh:
            item = curline.rstrip().split("\t")
            ec_counts.append(int(item[-1]))
            for ecitem in item[1:-1]:
                tidx, hidx, _ = tidx2coord[int(ecitem)]
                rowid[hidx].append(ecidx)
                colid[hidx].append(tidx)
            ecidx += 1

        for h in range(alnmat.num_haplotypes):
            alnmat.data.append(
                coo_matrix(
                    (np.ones(len(rowid[h])), (rowid[h], colid[h])),
                    shape=(alnmat.num_reads, alnmat.num_loci),
                    dtype=int,
                )
            )
        for h in range(alnmat.num_haplotypes):
            alnmat.data[h] = alnmat.data[h].tocsr()

        # alnmat = data[0]
        # for h in range(1, alnmat.num_haplotypes):
        #    alnmat = alnmat + ((2 ** h) * data[h])

        LOG.info("Creating EC count matrix")
        ec_counts = np.array(ec_counts)
        alnmat.count = csr_matrix(
            ec_counts
        ).transpose()  # Note: This is eventually a csc_matrix

        LOG.info("Creating objects of meta data")
        alnmat.lengths = np.zeros((alnmat.num_loci, alnmat.num_haplotypes), dtype=int)
        for _, tinfo in tidx2coord.items():
            alnmat.lengths[tinfo[0], tinfo[1]] = tinfo[2]
        # transcript_lengths[transcript_lengths < 1] = 1.0
        alnmat.hname = np.array(list(haplotype_idx.keys()))
        alnmat.lname = np.array(list(transcript_idx.keys()))

    LOG.info("Parsing complete")
    return alnmat


def convert(salmon_dir, ec_filename, sample, target_filename=None):
    LOG.debug("-------------------------------------------")
    LOG.debug("Parameters:")
    LOG.debug(f"  SALMON directory: {salmon_dir}")
    LOG.debug(f"  EC file: {ec_filename}")
    LOG.debug(f"  Sample: {sample}")
    LOG.debug(f"  Target file: {target_filename}")
    LOG.debug("-------------------------------------------")

    time0 = time.time()
    transcripts, haplotypes, alnmat, cntmat, transcript_lengths = parse_salmon_ec(
        salmon_dir, target_filename
    )
    LOG.info(
        f"{ec_filename} parsed in {utils.format_time(time0, time.time())}"
    )

    time1 = time.time()
    LOG.info(f"Converting and storing to {ec_filename}")
    bin_utils.ecsave(
        ec_filename,
        [sample],
        haplotypes,
        transcripts,
        transcript_lengths,
        alnmat,
        cntmat,
    )
    LOG.info(
        f"{ec_filename} created in {utils.format_time(time1, time.time())}"
    )


def convert2(salmon_dir, ec_filename, sample, target_filename=None):
    LOG.debug("-------------------------------------------")
    LOG.debug("Parameters:")
    LOG.debug(f"  SALMON directory: {salmon_dir}")
    LOG.debug(f"  EC file: {ec_filename}")
    LOG.debug(f"  Sample: {sample}")
    LOG.debug(f"  Target file: {target_filename}")
    LOG.debug("-------------------------------------------")

    time0 = time.time()
    # transcripts, haplotypes, alnmat, cntmat, transcript_lengths = parse_salmon_ec(salmon_dir, target_filename)
    alnmat = parse_salmon_ec(salmon_dir, target_filename)
    alnmat.sname = [sample]
    LOG.info(
        "{} parsed in {}".format(ec_filename, utils.format_time(time0, time.time()))
    )

    time1 = time.time()
    LOG.info(f"Converting and storing to {ec_filename}")
    # bin_utils.ecsave(ec_filename, [sample], haplotypes, transcripts, transcript_lengths, alnmat, cntmat)
    bin_utils.ecsave2(ec_filename, alnmat)
    LOG.info(
        "{} created in {}".format(ec_filename, utils.format_time(time1, time.time()))
    )
