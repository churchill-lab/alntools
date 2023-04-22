from collections import namedtuple
from collections import OrderedDict
from contextlib import closing
import multiprocessing
import os
import sqlite3
import struct
import time

from Bio import bgzf
import pysam

from alntools import utils

lock = None

LOG = utils.get_logger()

parse_fields = [
    "header_size",
    "begin_read_offset",
    "begin_read_size",
    "file_offset",
    "file_bytes",
    "end_read_offset",
    "end_read_size",
]
ParseRecord = namedtuple("ParseRecord", parse_fields)


class ParseParams(object):
    """
    # each core needs
    # - name of alignment file
    # - header size
    # - target file
    # - list of files to create and work on
    #   - idx, vo_start, vo_end

    """

    slots = ["input_file", "target_file", "temp_dir", "process_id", "data"]

    def __init__(self):
        self.input_file = None
        self.sql_file = None
        self.lock = None
        self.process_id = None
        self.data = []  # tuple of (idx, ParseRecord)

    def __str__(self):
        return (
            f"Input: {self.input_file}\n"
            f"Process ID: {self.process_id}\n"
            f"Data: {self.data}"
        )


class ParseResults(object):
    slots = ["read_ids"]

    def __init__(self):
        self.read_ids = None


def process_parse_fastq(cp):
    """

    :return:
    """
    LOG.debug(f"Process ID: {cp.process_id}, Input File: {cp.input_file}")

    temp_name = os.path.join(cp.temp_dir, "_fastq_parser.")

    rec_count = 0

    data = []

    SQL = """
    CREATE TABLE IF NOT EXISTS mapping (
       read_id TEXT,
       sample_id TEXT
    );
    """

    db_name = f"process_{cp.process_id}.dbs"

    try:
        for file_info_data in cp.data:
            reference_id = None
            reference_ids = []

            try:
                idx = file_info_data[0]
                parse_record = file_info_data[1]

                # must create the file
                temp_file = f"{temp_name}{idx}.bam"
                LOG.debug(
                    f"Process ID: {cp.process_id}, Creating alignment "
                    f"file: {temp_file}"
                )
                utils.delete_file(temp_file)
                chunk_fastq_file(cp.input_file, temp_file, parse_record)
                LOG.debug(
                    f"Process ID: {cp.process_id}, Opening alignment "
                    f"file: {temp_file}"
                )

                with pysam.FastxFile(temp_file) as fh:
                    rec_count = 0
                    for entry in fh:
                        rec_count += 1

                        data.append((entry.name, entry.sequence))
                        if rec_count % 100000 == 0:
                            # lock.acquire()
                            with closing(sqlite3.connect(db_name)) as conn:
                                cursor = conn.cursor()
                                cursor.execute(SQL)
                                cursor.executemany(
                                    "insert into mapping values (?,?)", data
                                )
                                conn.commit()
                            # lock.release()

                            data = []
                            LOG.debug(
                                f"Process ID: {cp.process_id}, File: "
                                f"{temp_file}, {rec_count:,} "
                                "Fastx records processed"
                            )

                # lock.acquire()
                with closing(sqlite3.connect(db_name)) as conn:
                    cursor = conn.cursor()
                    cursor.execute(SQL)
                    cursor.executemany("insert into mapping values (?,?)", data)
                    conn.commit()
                    data = []
                # lock.release()

            except Exception as e1:
                print("ERROR")
                print(str(e1))

            LOG.info(
                f"DONE Process ID: {cp.process_id}, File: {temp_file}, "
                f"{rec_count:,} Fastx records processed"
            )

            utils.delete_file(temp_file)

    except Exception as e:
        LOG.error(f"Error: {e}")

    ret = ParseResults()
    ret.read_ids = []

    return ret


def wrapper_convert(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    # print str(args)
    return process_parse_fastq(*args)


def init(l):
    global lock
    lock = l


# barcode_file
# read_id
# sequence
# +
# quality
#
# read_id line = @ST-K00126:296:HCK57BBXX:7:1101:2889:1156 1:N:0:ACAGAGGT
# read_id = ST-K00126:296:HCK57BBXX:7:1101:2889:1156


def parse(fastq_filename, num_chunks=0, temp_dir=None):
    """ """
    start_time = time.time()

    num_processes = multiprocessing.cpu_count()

    if num_chunks <= 0:
        num_chunks = num_processes
    else:
        if num_chunks > 1000:
            LOG.info(f"Modifying number of chunks from {num_chunks} to 1000")
            num_chunks = 1000

    if not temp_dir:
        temp_dir = os.path.dirname(fastq_filename)

    LOG.info(f"Calculating {num_chunks:,} chunks")
    temp_time = time.time()
    chunks = calculate_chunks(fastq_filename, num_chunks)
    LOG.info(
        f"{len(chunks):,} chunks calculated in "
        f"{utils.format_time(temp_time, time.time())}, "
        f"total time: {utils.format_time(start_time, time.time())}"
    )

    # each core needs

    l = multiprocessing.Lock()

    all_params = []

    pid = 0
    for temp_chunk_ids in utils.partition(
        [idx for idx in range(num_chunks)], num_processes
    ):
        params = ParseParams()
        params.input_file = fastq_filename
        params.temp_dir = temp_dir

        for x, cid in enumerate(temp_chunk_ids):
            params.process_id = pid
            params.data.append((cid, chunks[cid]))

        pid += 1
        all_params.append(params)
        LOG.debug(f"params = {params}")

    final = ParseResults()
    final.read_ids = OrderedDict()

    LOG.info(f"Starting {num_processes} processes ...")

    temp_time = time.time()
    args = zip(all_params)
    pool = multiprocessing.Pool(
        initializer=init, initargs=(l,), processes=num_processes
    )
    results = pool.map(wrapper_convert, args)

    LOG.info(
        f"All processes done in {utils.format_time(temp_time, time.time())}, "
        f"total time: {utils.format_time(start_time, time.time())}"
    )


def split_fastq(filename, n, directory=None):
    """
    Split a FASTQ file into ``n`` files.

    :param str filename: the name of the FASTQ file
    :param int n: number of files to chunk into
    :param str directory: output directory, defaults to ``bam_filename`` directory

    :return: names of the files
    """
    start_time = time.time()

    LOG.debug(f"FASTQ File: {filename}")
    LOG.debug(f"Number of Files: {n}")

    if not directory:
        directory = os.path.dirname(filename)

    LOG.debug(f"Output Directory: {directory}")

    bam_basename = os.path.basename(filename)
    bam_prefixname, bam_extension = os.path.splitext(bam_basename)
    bam_output_temp = os.path.join(directory, bam_prefixname)

    LOG.info(f"Calculating {n:,} chunks...")
    temp_time = time.time()
    chunks = calculate_chunks(filename, n)
    LOG.info(
        f"{len(chunks):,} chunks calculated in "
        f"{utils.format_time(temp_time, time.time())}, "
        f"total time: {utils.format_time(start_time, time.time())}"
    )

    file_names = []

    for idx, chunk in enumerate(chunks):
        # must create the file
        LOG.debug(chunk)
        new_file = f"{bam_output_temp}_{idx}{bam_extension}"
        file_names.append(new_file)
        LOG.debug(f"Creating FASTQ file: {new_file}")
        chunk_fastq_file(filename, new_file, chunk)

    LOG.info(
        f"{len(chunks):,} files created in "
        f"{utils.format_time(temp_time, time.time())}, "
        f"total time: {utils.format_time(start_time, time.time())}"
    )

    return file_names


def bytes_from_file(read_filename, write_filename, offset=0, bytes_size=-1):
    """
    Read bytes from a file and append them onto another file.

    :param read_filename: the name of the file to read from
    :param write_filename: the name of the file to write to
    :param offset: the number of bytes to offset (seek)
    :param bytes_size: the number of bytes to read, -1 = to end of file
    :return:
    """
    with open(read_filename, "rb") as fr:
        if offset > 0:
            fr.seek(offset)

        if bytes_size > 0:
            # specific size
            with open(write_filename, "ab") as fw:
                fw.write(fr.read(bytes_size))
        elif bytes_size < 0:
            # entire file
            with open(write_filename, "ab") as fw:
                fw.write(fr.read())
        else:
            # just creating an empty file
            open(write_filename, "a").close()


def chunk_fastq_file(fastq_filename, new_filename, parse_rec):
    """
    Create a new FASTQ file from an existing one.

    :param str fastq_filename: the name of the original BAM file
    :param str new_filename: the name of the new BAM file
    :param class:`ParseRecord` parse_rec: the information containing where to extract
    :return:
    """
    try:
        os.remove(new_filename)
    except Exception as e:
        pass

    # copy the header from original BAM file to new
    bytes_from_file(fastq_filename, new_filename, 0, parse_rec.header_size)

    if parse_rec.begin_read_offset > 0:
        # if there are reads before a chunk offset, we need to extract them
        b = bgzf.BgzfReader(fastq_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.begin_read_offset)
        b2.write(b.read(parse_rec.begin_read_size))
        b2.close()

    # grab bgzf chunks from the OLD FASTQ file and append to NEW FASTQ file
    bytes_from_file(
        fastq_filename, new_filename, parse_rec.file_offset, parse_rec.file_bytes
    )

    if parse_rec.end_read_offset > 0:
        # if there are reads after a chunk offset, we need to extract them
        b = bgzf.BgzfReader(fastq_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.end_read_offset)
        b2.write(b.read(parse_rec.end_read_size))
        b2.close()


def calculate_chunks(filename, num_chunks):
    """
    Calculate the boundaries in the BAM file and partition into chunks.

    :param str filename: name of the BAM file
    :param int num_chunks: number of chunks to partition the boundaries into
    :return: a list of tuples containing the start and end boundaries
    """
    if num_chunks <= 0:
        raise ValueError("The number of chunks to calculate should be >= 1")

    if num_chunks == 1:
        # aln_file = pysam.AlignmentFile(filename)
        # header_size = bgzf.split_virtual_offset(aln_file.tell())[0]
        # aln_file.close()

        pr = ParseRecord(0, 0, 0, 0, -1, 0, 0)
        return [pr]

    try:
        f = open(filename, "r")
        # get all the block start offsets
        block_offsets = []
        decompressed_lengths = []
        i = 0

        for values in FastBgzfBlocks(f):
            block_offsets.append(values[0])
            decompressed_lengths.append(values[3])
            i += 1

        # partition the starts into manageable chunks
        div, mod = divmod(len(block_offsets), num_chunks)

        fastq_fh = bgzf.BgzfReader(filename, "r")
        header_size = 0
        partitioned_offsets = [(header_size, 0)]

        for i in range(1, num_chunks):
            index = div * i + min(i, mod)

            virtual_offset = bgzf.make_virtual_offset(block_offsets[index], 0)
            fastq_fh.seek(virtual_offset)
            line = fastq_fh.readline().strip()
            while line != "+":
                line = fastq_fh.readline().strip()
            quality_line = fastq_fh.readline()
            virtual_offset = fastq_fh.tell()

            # block start & within block offset
            partitioned_offsets.append(bgzf.split_virtual_offset(virtual_offset))

        fastq_fh.close()

        # now let's calculate beginning and ends
        params = []

        for i, offset in enumerate(partitioned_offsets):
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
                # print 'file_offset=', file_offset
                file_bytes = partitioned_offsets[i + 1][0] - file_offset
                # print 'file_bytes=', file_bytes
                end_read_offset = bgzf.make_virtual_offset(
                    partitioned_offsets[i + 1][0], 0
                )
                end_read_size = partitioned_offsets[i + 1][1]
            elif i == num_chunks - 1:
                # last
                begin_read_offset = bgzf.make_virtual_offset(
                    partitioned_offsets[i][0], partitioned_offsets[i][1]
                )
                begin_read_size = (
                    decompressed_lengths[index] - partitioned_offsets[i][1]
                )
                file_offset = block_offsets[index + 1]
                file_bytes = -1
                end_read_offset = 0
                end_read_size = 0
            else:
                # all others
                if offset[1] == 0:
                    # bgzf boundary
                    print("****************HUH")
                    return

                begin_read_offset = bgzf.make_virtual_offset(
                    partitioned_offsets[i][0], partitioned_offsets[i][1]
                )
                begin_read_size = (
                    decompressed_lengths[index] - partitioned_offsets[i][1]
                )
                file_offset = block_offsets[index + 1]
                file_bytes = partitioned_offsets[i + 1][0] - file_offset

                end_read_offset = bgzf.make_virtual_offset(
                    partitioned_offsets[i + 1][0], 0
                )
                end_read_size = partitioned_offsets[i + 1][1]

            pr = ParseRecord(
                header_size,
                begin_read_offset,
                begin_read_size,
                file_offset,
                file_bytes,
                end_read_offset,
                end_read_size,
            )
            params.append(pr)

        return params

    except Exception as e:
        print(f"calculate_chunks error: {e}")


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
        raise ValueError(
            "A BGZF block should start with %r, not %r; handle.tell() now says %r"
            % (bgzf._bgzf_magic, magic, handle.tell())
        )

    gzip_mod_time, gzip_extra_flags, gzip_os, extra_len = struct.unpack(
        "<LBBH", handle.read(8)
    )
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


if __name__ == "__main__":
    parse("data/9000samp.readID2sampID.fastq.gz", 10, "tmp")
    # 'CREATE INDEX IF NOT EXISTS idx_mapping_read_id ON mapping (read_id ASC)',
    # 'CREATE INDEX IF NOT EXISTS idx_mapping_sample_id ON mapping (sample_id ASC)',
