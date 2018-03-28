from collections import namedtuple
from Bio import bgzf
import os


from alntools import utils
from alntools.bam_utils import FastBgzfBlocks

parse_fields = ["header_size", "begin_read_offset", "begin_read_size",
                "file_offset", "file_bytes", "end_read_offset",
                "end_read_size"]

ParseRecord = namedtuple("ParseRecord", parse_fields)

import multiprocessing

def truncate_fastq_file(fname):
    """
    Remove the EOF from BGZF/BAM file.

    Does not check if the EOF is present or not.

    :param fname: the name of the BZF/BAM file
    :return:
    """
    utils.truncate_file(fname, 28)



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
            open(write_filename, 'a').close()



def chunk_fastq_file(fastq_filename, new_filename, parse_rec):
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
    bytes_from_file(fastq_filename, new_filename, 0, parse_rec.header_size)

    if parse_rec.begin_read_offset > 0:
        # if there are reads before a chunk offset, we need to extract them
        b = bgzf.BgzfReader(fastq_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.begin_read_offset)
        b2.write(b.read(parse_rec.begin_read_size))
        b2.close()
        #truncate_fastq_file(new_filename)

    # grab bgzf chunks from the OLD BAM file and append to NEW BAM file
    bytes_from_file(fastq_filename, new_filename, parse_rec.file_offset, parse_rec.file_bytes)

    if parse_rec.end_read_offset > 0:
        # if there are reads after a chunk offset, we need to extract them
        b = bgzf.BgzfReader(fastq_filename)
        b2 = bgzf.BgzfWriter(new_filename, mode="a")
        b.seek(parse_rec.end_read_offset)
        b2.write(b.read(parse_rec.end_read_size))
        b2.close()

    # fix the bam EOF if needed
    #fix_bam(new_filename)


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
        f = open(filename, 'r')
        # get all the block start offsets
        block_offsets = []
        decompressed_lengths = []
        i = 0

        for values in FastBgzfBlocks(f):
            block_offsets.append(values[0])
            decompressed_lengths.append(values[3])

            if i % 50000 == 0:
                print  'Block {}'.format(i)
            i += 1

        # partition the starts into manageable chunks
        div, mod = divmod(len(block_offsets), num_chunks)

        fastq_fh = bgzf.BgzfReader(filename, 'r')
        header_size = 0
        partitioned_offsets = [(header_size, 0)]

        for i in range(1, num_chunks):
            index = div * i + min(i, mod)

            virtual_offset = bgzf.make_virtual_offset(block_offsets[index], 0)
            print 'block offsets = ', block_offsets[index]
            print 'virtual_offset = ', virtual_offset
            fastq_fh.seek(virtual_offset)
            line = fastq_fh.readline().strip()
            while line != '+':
                line = fastq_fh.readline().strip()
            quality_line = fastq_fh.readline()
            print 'quality_line=', quality_line
            virtual_offset = fastq_fh.tell()

            # block start & within block offset
            partitioned_offsets.append(
                bgzf.split_virtual_offset(virtual_offset))

        fastq_fh.close()

        print 'partitioned_offsets=', str(partitioned_offsets)


        # now let's calculate beginning and ends
        params = []

        for i, offset in enumerate(partitioned_offsets):
            print '{} => {}'.format(i, offset)

            index = block_offsets.index(partitioned_offsets[i][0])
            print 'index=', str(index)
            begin_read_offset = 0
            begin_read_size = 0
            file_offset = 0
            file_bytes = 0
            end_read_offset = 0
            end_read_size = 0

            if i == 0:
                # first
                print 'first'
                begin_read_offset = 0
                begin_read_size = 0
                file_offset = block_offsets[index]
                # print 'file_offset=', file_offset
                file_bytes = partitioned_offsets[i + 1][0] - file_offset
                # print 'file_bytes=', file_bytes
                end_read_offset = bgzf.make_virtual_offset(
                        partitioned_offsets[i + 1][0], 0)
                end_read_size = partitioned_offsets[i + 1][1]
            elif i == num_chunks - 1:
                # last
                print 'last'
                begin_read_offset = bgzf.make_virtual_offset(
                        partitioned_offsets[i][0], partitioned_offsets[i][1])
                begin_read_size = decompressed_lengths[index] - \
                                  partitioned_offsets[i][1]
                file_offset = block_offsets[index + 1]
                file_bytes = -1
                end_read_offset = 0
                end_read_size = 0
            else:
                # all others
                print 'other'
                if offset[1] == 0:
                    # bgzf boundary
                    print '****************HUH'
                    return

                print 'calculating begin read offset from partitioned_offsets[i] = ', str(partitioned_offsets[i])
                begin_read_offset = bgzf.make_virtual_offset(
                        partitioned_offsets[i][0], partitioned_offsets[i][1])

                print 'calculating begin read size from decompressed_lengths[index] and partioned offsets',  decompressed_lengths[index], str(partitioned_offsets[i][1])

                begin_read_size = decompressed_lengths[index] - \
                                  partitioned_offsets[i][1]

                print 'calculating file offset from block_offsets[index+1] = ', block_offsets[index + 1]

                file_offset = block_offsets[index + 1]

                print 'calculating file_bytes from partitioned_offsets[[i + 1][0] - file_offset= ',  partitioned_offsets[i + 1][0], file_offset, (partitioned_offsets[i + 1][0] - file_offset)

                file_bytes = partitioned_offsets[i + 1][0] - file_offset


                end_read_offset = bgzf.make_virtual_offset(
                        partitioned_offsets[i + 1][0], 0)
                end_read_size = partitioned_offsets[i + 1][1]

            pr = ParseRecord(header_size, begin_read_offset, begin_read_size,
                             file_offset, file_bytes, end_read_offset,
                             end_read_size)
            params.append(pr)

        return params

    except Exception as e:
        print 'calculate_chunks error: {}'.format(str(e))


def make_more_fastq(file_name, num):
    chunks = calculate_chunks(file_name, num)



    for i, c in enumerate(chunks):
        print str(i), str(c)



    for i, c in enumerate(chunks):
        print 'making ' +  '{}.vcf'.format(i)
        chunk_fastq_file(file_name, '{}.vcf'.format(i), c)


if __name__ == '__main__':
    file_name = "data/9000samp.readID2sampID.fastq.gz"
    make_more_fastq(file_name, 4)
