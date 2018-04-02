# -*- coding: utf-8 -*-
#
# TODO: CHANGE LOGGING
#
# I DON'T LIKE THIS METHOD OF LOGGING
#
# I WANTED A VERBOSE METHOD, BUT DIDN'T WANT TO ADD A NEW LEVEL
#
# THUS...
#
# logging.WARNING is informational
# logging.INFO is user debug
# logging.DEBUG is developer debug
#
from collections import OrderedDict
from past.builtins import xrange

import logging
import os

logging.basicConfig(format='[alntools] [%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


def get_logger():
    """Get the :class:`logging.Logger`.

    Returns:
        :class:`logging.Logger`
    """
    return logging.getLogger(__name__)


def configure_logging(level):
    """Configure the :class:`Logger`.

    - 0 = WARNING
    - 1 = INFO
    - 2 = DEBUG

    Args:
        level (int): The level to debug.
    """
    if level == 0:
        get_logger().setLevel(logging.WARN)
    elif level == 1:
        get_logger().setLevel(logging.INFO)
    elif level > 1:
        get_logger().setLevel(logging.DEBUG)


def format_time(start, end):
    """Format length of time between ``start`` and ``end``.

    Args:
        start (:class:`time.time`): The start time.
        end (:class:`time.time`): The end time.

    Returns:
        str: A formatted string of hours, minutes, and seconds.

    """
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def partition(lst, n):
    """Split ``lst`` into ``n`` partitions.

    Args:
        lst (list): A list of values to be partitioned.
        n (int): The number of partitions.

    Returns:
        list: A list of length n, each containing a sublist.
    """
    q, r = divmod(len(lst), n)
    indices = [q*i + min(i, r) for i in xrange(n+1)]
    temp = [lst[indices[i]:indices[i+1]] for i in xrange(n)]

    ret = []
    for x in temp:
        if len(x) == 0:
            break
        ret.append(x)

    return ret


def list_to_int(lst):
    """Convert a list (assuming values of 0 or 1) into an integer

    Example: [0, 0, 0] => 0  (1 * 0 + 2 * 0 + 3 * 0)
             [0, 1, 0] => 2  (1 * 0 + 2 * 1 + 3 * 0)
             [1, 1, 0] => 3  (1 * 1 + 2 * 1 + 3 * 0)
             [1, 1, 1] => 7  (1 * 1 + 2 * 1 + 3 * 1)
             [1, 0, 0, 1, 1, 1, 1] => 121

    Args:
        lst (list): The list.

    Returns:
        int: The decimal representation.
    """
    c = 0
    for i, on_off in enumerate(lst):
        if on_off == 1:
            c |= 1 << i
    return c


def int_to_list(c, size):
    """Convert an integer to a list on length size

    Example: [0, 0, 0] => 0  (1 * 0 + 2 * 0 + 3 * 0)
             [0, 1, 0] => 2  (1 * 0 + 2 * 1 + 3 * 0)
             [1, 1, 0] => 3  (1 * 1 + 2 * 1 + 3 * 0)
             [1, 1, 1] => 7  (1 * 1 + 2 * 1 + 3 * 1)
             [1, 0, 0, 1, 1, 1, 1] => 121

    Args:
        c (int): The integer value.
        size: The size of the list to be returned.

    Returns:
        list: A list of size ``n`` representing the value ``c``.

    """
    ret = [0] * size
    for i in xrange(0, size):
        if (c & (1 << i)) != 0:
            ret[i] = 1
    return ret


def bytes_from_file(read_filename, write_filename, offset=0, bytes_size=-1):
    """Read bytes from a file and append them onto another file.

    Args:
        read_filename (str): The name of the file to read from.
        write_filename (str): The name of the file to write to.
        offset (int): The number of bytes to offset (seek).
        bytes_size: The number of bytes to read, -1 = to end of file.

    Returns:
        None
    """
    with open(read_filename, "rb") as fr:
        if offset > 0:
            fr.seek(offset)

        if bytes_size > 0:
            data = fr.read(bytes_size)
        else:
            data = fr.read()

        with open(write_filename, "ab") as fw:
            fw.write(data)


def parse_targets(target_file):
    """Parse a file and return the targets.

    Args:
        target_file (str): The name of the file.

    Returns:
        OrderedDict: The target name as the key and the insertion order as the
            value.
    """
    targets = OrderedDict()
    with open(target_file, 'r') as f:
        for line in f:
            if line and line[0] == '#':
                continue
            _id = line.strip().split()[0]
            targets[_id] = len(targets)
    return targets


def delete_file(file_name):
    """Delete a file.

    Args:
        file_name (str): The file name.
    """
    try:
        os.remove(file_name)
    except:
        pass


def truncate_file(file_name, bytes_from_end):
    """Truncate a file.

    Args:
        file_name (str): The name of the file.
        bytes_from_end (int): Number of bytes from end of file.
    """
    f = open(file_name, 'r+')
    f.seek(-1 * bytes_from_end, os.SEEK_END)
    f.truncate()
    f.close()


def get_bam_files(files):
    """Get all BAM files in tuple of files

    Args:
        files (tuple): A tuple of files or directories.

    Returns:
        list: A list of files in directories.
    """
    list_files = []

    for file in files:
        if os.path.isfile(file):
            list_files.append(file)
        elif os.path.isdir(file):
            list_files.extend(get_files_in_dir(file, ['bam']))

    return list_files


def get_files_in_dir(directory, file_extensions=None):
    """Get a list of all files in ``directory``.

    Args:
        directory (str): The directory to search.
        file_extensions: A list of file extensions to match, None means all.

    Returns:
        list: A list of all matching files.
    """
    matches = []

    for root, directories, filenames in os.walk(directory):
        for filename in filenames:
            if file_extensions:
                for ext in file_extensions:
                    if filename.lower().endswith(ext.lower()):
                        matches.append(os.path.join(root, filename))
            else:
                matches.append(os.path.join(root, filename))
        # for directory in directories:
        #    print(os.path.join(root, directory))

    return matches


