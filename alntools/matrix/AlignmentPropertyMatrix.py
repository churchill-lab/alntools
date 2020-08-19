#!/usr/bin/env python
import copy
import tables
from struct import pack, unpack
from itertools import dropwhile

import numpy as np
from scipy.sparse import lil_matrix, coo_matrix, csc_matrix, csr_matrix
from alntools.matrix.Sparse3DMatrix import Sparse3DMatrix

try:
    xrange
except NameError:
    xrange = range


def enum(**enums):
    return type('Enum', (), enums)


def is_comment(s):
    return s.startswith('#')


class AlignmentPropertyMatrix(Sparse3DMatrix):

    Axis = enum(LOCUS=0, HAPLOTYPE=1, READ=2, GROUP=3, HAPLOGROUP=4)

    def __init__(self, 
                 other=None,
                 ecfile=None,
                 h5file=None, datanode='/', metanode='/', shallow=False,
                 shape=None, dtype=float, haplotype_names=None, locus_names=None, read_names=None, sample_names=None,
                 grpfile=None):

        Sparse3DMatrix.__init__(self, other=other, h5file=h5file, datanode=datanode, shape=shape, dtype=dtype)

        self.num_loci, self.num_haplotypes, self.num_reads = self.shape
        self.num_samples = 0
        self.num_groups = 0
        self.count   = None
        self.hname   = None
        self.lname   = None   # locus name
        self.rname   = None   # read name
        self.sname   = None   # sample name (e.g. sample barcodes::cell barcodes)
        self.lid     = None   # locus ID
        self.rid     = None   # read ID
        self.sid     = None   # sample ID
        self.lengths = None   # transcript lengths (or effective lengths)
        self.gname   = None   # group name
        self.groups  = None   # groups in terms of locus IDs

        if other is not None:  # Use for copying from other existing AlignmentPropertyMatrix object
            if other.count is not None:
                self.count = copy.copy(other.count)
            if other.lengths is not None:
                self.lengths = copy.copy(other.lengths)
            if not shallow:
                self.__copy_names(other)
                self.__copy_group_info(other)

        elif ecfile is not None:
            with open(ecfile, 'rb') as f:
                ecformat = unpack('<i', f.read(4))[0]
                if ecformat == 2:
                    self.num_haplotypes = unpack('<i', f.read(4))[0]
                    hname = list()
                    for hidx in range(self.num_haplotypes):
                        hname_len = unpack('<i', f.read(4))[0]
                        hname.append(unpack('<{}s'.format(hname_len), f.read(hname_len))[0].decode('utf-8'))
                    self.hname = np.array(hname)
                    # Get transcipt info
                    self.num_loci = unpack('<i', f.read(4))[0]
                    self.lengths = np.zeros((self.num_loci, self.num_haplotypes), dtype=float)
                    tname = list()
                    for tidx in range(self.num_loci):
                        tname_len = unpack('<i', f.read(4))[0]
                        tname.append(unpack('<{}s'.format(tname_len), f.read(tname_len))[0].decode('utf-8'))
                        for hidx in range(self.num_haplotypes):
                            self.lengths[tidx, hidx] = unpack('<i', f.read(4))[0]
                    self.lname = np.array(tname)
                    self.lid = dict(zip(self.lname, np.arange(self.num_loci)))
                    # Get sample info
                    sname = list()
                    self.num_samples = unpack('<i', f.read(4))[0]
                    for sidx in range(self.num_samples):
                        sname_len = unpack('<i', f.read(4))[0]
                        sname.append(unpack('<{}s'.format(sname_len), f.read(sname_len))[0].decode('utf-8'))
                    self.sname = np.array(sname)
                    self.sid = dict(zip(self.sname, np.arange(self.num_samples)))
                    # Read in alignment matrix info
                    indptr_len = unpack('<i', f.read(4))[0]
                    self.num_reads = indptr_len - 1
                    nnz = unpack('<i', f.read(4))[0]
                    indptr_A = np.array(unpack('<{}i'.format(indptr_len), f.read(4*indptr_len)))
                    indices_A = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
                    data_A = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
                    # Read in EC count matrix info
                    indptr_len = unpack('<i', f.read(4))[0]
                    nnz = unpack('<i', f.read(4))[0]
                    indptr_N = np.array(unpack('<{}i'.format(indptr_len), f.read(4*indptr_len)))
                    indices_N = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
                    data_N = np.array(unpack('<{}i'.format(nnz), f.read(4*nnz)))
                    # Populate class member variables
                    for hidx in range(self.num_haplotypes-1):
                        data_A, data_A_rem = np.divmod(data_A, 2)
                        self.data.append(csr_matrix((data_A_rem, indices_A, indptr_A), shape=(self.num_reads, self.num_loci)))
                    self.data.append(csr_matrix((data_A, indices_A, indptr_A), shape=(self.num_reads, self.num_loci)))
                    for hidx in range(self.num_haplotypes):
                        self.data[hidx].eliminate_zeros()
                    self.count = csc_matrix((data_N, indices_N, indptr_N), shape=(self.num_reads, self.num_samples))
                    if self.num_samples == 1:
                        self.count = self.count.todense().A.flatten()
                    self.shape = (self.num_loci, self.num_haplotypes, self.num_reads)
                    self.finalize()
                elif ecformat == 1:
                    raise NotImplementedError
                elif ecformat == 0:
                    raise TypeError('Format 0 is not supported anymore.')
 
        elif h5file is not None:  # Use for loading from a pytables file
            h5fh = tables.open_file(h5file, 'r')
            if not shallow:
                self.hname = h5fh.get_node_attr(datanode, 'hname')
                self.lname = np.char.decode(h5fh.get_node(metanode, 'lname').read(), 'utf-8')
                self.lid = dict(zip(self.lname, np.arange(self.num_loci)))
                if h5fh.__contains__('%s' % (metanode + '/rname')):
                    self.rname = np.char.decode(h5fh.get_node(metanode, 'rname').read(), 'utf-8')
                    self.rid = dict(zip(self.rname, np.arange(self.num_reads)))
                if h5fh.__contains__('%s' % (metanode + '/sname')):
                    self.sname = np.char.decode(h5fh.get_node(metanode, 'sname').read(), 'utf-8')
                    self.sid = dict(zip(self.sname, np.arange(self.num_samples)))
                    self.num_samples = len(self.sname)
            if h5fh.__contains__('%s' % (datanode + '/count')):
                try:
                    self.count = h5fh.get_node(datanode, 'count').read()  # Format-1
                except tables.NoSuchNodeError as e:  # Format-2
                    nmat_node = h5fh.get_node(datanode + '/count')
                    indptr = h5fh.get_node(nmat_node, 'indptr').read()
                    indices = h5fh.get_node(nmat_node, 'indices').read()
                    data = h5fh.get_node(nmat_node, 'data').read()
                    self.count = csc_matrix((data, indices, indptr), dtype=np.float64)
                    self.num_samples = self.count.shape[1]
            if h5fh.__contains__('%s' % (datanode + '/lengths')):
                self.lengths = h5fh.get_node(datanode, 'lengths').read()
            h5fh.close()

        elif shape is not None:  # Use for initializing an empty matrix
            if haplotype_names is not None:
                if len(haplotype_names) == self.num_haplotypes:
                    self.hname = haplotype_names
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')
            if locus_names is not None:
                if len(locus_names) == self.num_loci:
                    self.lname = np.array(locus_names)
                    self.lid   = dict(zip(self.lname, np.arange(self.num_loci)))
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')
            if read_names is not None:
                if len(read_names) == self.num_reads:
                    self.rname = np.array(read_names)
                    self.rid   = dict(zip(self.rname, np.arange(self.num_reads)))
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')
            if sample_names is not None:
                self.sname = np.array(sample_names)
                self.sid   = dict(zip(self.sname, np.arange(self.num_samples)))
                self.num_samples = len(sample_names)
            else:
                self.num_samples = 1

        if grpfile is not None:
            self.__load_groups(grpfile)

    def __load_groups(self, grpfile):  # A group is a set of isoforms within a gene
        if self.lid is not None:
            self.gname  = list()
            self.groups = list()
            with open(grpfile) as fh:
                for curline in fh:
                    item = curline.rstrip().split("\t")
                    self.gname.append(item[0])
                    tid_list = [ self.lid[t] for t in item[1:] ]
                    self.groups.append(tid_list)
            self.gname = np.array(self.gname)
            self.num_groups = len(self.gname)
        else:
            raise RuntimeError('Locus IDs are not availalbe.')

    load_groups = __load_groups

    def __copy_names(self, other):
        self.hname = other.hname
        self.lname = copy.copy(other.lname)
        self.rname = copy.copy(other.rname)
        self.sname = copy.copy(other.sname)
        self.lid   = copy.copy(other.lid)
        self.rid   = copy.copy(other.rid)
        self.sid   = copy.copy(other.sid)

    def __copy_group_info(self, other):
        # if other.num_groups > 0:
        if other.groups is not None and other.gname is not None:
            self.groups = copy.deepcopy(other.groups)
            self.gname  = copy.copy(other.gname)
            self.num_groups = other.num_groups

    def copy(self, shallow=False):
        dmat = Sparse3DMatrix.copy(self)
        dmat.count = self.count.copy()
        dmat.lengths = self.lengths.copy()
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape
        if not shallow:
            dmat.__copy_names(self)
            dmat.__copy_group_info(self)
        return dmat

    def _bundle_inline(self, reset=False): # Inline bundling method
        if self.finalized:
            if self.num_groups > 0 and self.groups is not None and self.gname is not None:
                grp_conv_mat = lil_matrix((self.num_loci, self.num_groups))
                for i in xrange(self.num_groups):
                    grp_conv_mat[self.groups[i], i] = 1.0
                grp_conv_mat = grp_conv_mat.tocsc()
                for hid in xrange(self.num_haplotypes):
                    self.data[hid] = self.data[hid] * grp_conv_mat # TODO: Is there any better way to save memory?
                self.num_loci = self.num_groups
                self.shape = (self.num_groups, self.num_haplotypes, self.num_reads)
                self.lname = copy.copy(self.gname)
                self.lid   = dict(zip(self.gname, np.arange(self.num_groups)))
                self.num_groups = 0
                self.groups = None
                self.gname = None
                if reset:
                    self.reset()
            else:
                raise RuntimeError('No group information is available for bundling.')
        else:
            raise RuntimeError('The matrix is not finalized.')

    def bundle(self, reset=False, shallow=False): # Copies the original matrix (Use lots of memory)
        """
        Returns ``AlignmentPropertyMatrix`` object in which loci are bundled using grouping information.

        :param reset: whether to reset the values at the loci
        :param shallow: whether to copy all the meta data
        """
        if self.finalized:
            # if self.num_groups > 0:
            if self.groups is not None and self.gname is not None:
                grp_conv_mat = lil_matrix((self.num_loci, self.num_groups))
                for i in xrange(self.num_groups):
                    grp_conv_mat[self.groups[i], i] = 1.0
                grp_align = Sparse3DMatrix.__mul__(self, grp_conv_mat) # The core of the bundling
                grp_align.num_loci = self.num_groups
                grp_align.num_haplotypes = self.num_haplotypes
                grp_align.num_reads = self.num_reads
                grp_align.shape = (grp_align.num_loci, grp_align.num_haplotypes, grp_align.num_reads)
                grp_align.count = self.count
                if not shallow:
                    grp_align.lname = copy.copy(self.gname)
                    grp_align.hname = self.hname
                    grp_align.rname = copy.copy(self.rname)
                    grp_align.sname = copy.copy(self.sname)
                    grp_align.lid   = dict(zip(grp_align.lname, np.arange(grp_align.num_loci)))
                    grp_align.rid   = copy.copy(self.rid)
                    grp_align.sid   = copy.copy(self.sid)
                if reset:
                    grp_align.reset()
                return grp_align
            else:
                raise RuntimeError('No group information is available for bundling.')
        else:
            raise RuntimeError('The matrix is not finalized.')

    #
    # Binary Operators
    #

    def __add__(self, other):
        dmat = Sparse3DMatrix.__add__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = self.shape
        dmat.__copy_names(self)
        dmat.__copy_group_info(self)
        return dmat

    def __sub__(self, other):
        dmat = Sparse3DMatrix.__sub__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = self.shape
        dmat.__copy_names(self)
        dmat.__copy_group_info(self)
        return dmat

    def __mul__(self, other):
        dmat = Sparse3DMatrix.__mul__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape
        if isinstance(other, (np.ndarray, csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
            dmat.hname = self.hname
            dmat.rname = copy.copy(self.rname)
            dmat.rid   = copy.copy(self.rid)
            dmat.num_groups = 0
        else:
            dmat.__copy_names(self)
            dmat.__copy_group_info(self)
        return dmat

    #
    # Helper functions
    #

    def sum(self, axis):
        if self.finalized:
            if axis == self.Axis.LOCUS:
                sum_mat = []  # sum along loci
                for hid in xrange(self.num_haplotypes):
                    sum_mat.append(self.data[hid].sum(axis=1).A)
                sum_mat = np.hstack(sum_mat)
            elif axis == self.Axis.HAPLOTYPE:  # sum along haplotypes
                sum_mat = self.data[0]
                for hid in xrange(1, self.num_haplotypes):
                    sum_mat = sum_mat + self.data[hid]  # Unlike others, this sum_mat is still sparse matrix
            elif axis == self.Axis.READ:  # sum along reads
                sum_mat = []
                for hid in xrange(self.num_haplotypes):
                    if self.count is None:
                        sum_hap = self.data[hid].sum(axis=0).A
                    else:
                        hap_mat = self.data[hid].copy()
                        hap_mat.data *= self.count[hap_mat.indices]
                        sum_hap = hap_mat.sum(axis=0).A
                    sum_mat.append(sum_hap)
                sum_mat = np.vstack(sum_mat)
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
            return sum_mat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def normalize_reads(self, axis, grouping_mat=None):
        """
        Read-wise normalization

        :param axis: The dimension along which we want to normalize values
        :param grouping_mat: An incidence matrix that specifies which isoforms are from a same gene
        :return: Nothing (as the method performs in-place operations)
        :rtype: None
        """
        if self.finalized:
            if axis == self.Axis.LOCUS:  # Locus-wise normalization on each read
                normalizer = self.sum(axis=self.Axis.HAPLOTYPE)  # Sparse matrix of |reads| x |loci|
                normalizer.eliminate_zeros()
                for hid in xrange(self.num_haplotypes):
                    self.data[hid].eliminate_zeros()  # Trying to avoid numerical problem (inf or nan)
                    self.data[hid] = np.divide(self.data[hid], normalizer)  # element-wise division
            elif axis == self.Axis.HAPLOTYPE:  # haplotype-wise normalization on each read
                for hid in xrange(self.num_haplotypes):
                    normalizer = self.data[hid].sum(axis=self.Axis.HAPLOTYPE)  # 1-dim Sparse matrix of |reads| x 1
                    normalizer = normalizer.A.flatten()
                    self.data[hid].data /= normalizer[self.data[hid].indices]
            elif axis == self.Axis.READ:  # normalization each read as a whole
                sum_mat = self.sum(axis=self.Axis.LOCUS)
                normalizer = sum_mat.sum(axis=self.Axis.HAPLOTYPE)
                normalizer = normalizer.ravel()
                for hid in xrange(self.num_haplotypes):
                    self.data[hid].data /= normalizer[self.data[hid].indices]
            elif axis == self.Axis.GROUP:  # group-wise normalization on each read
                if grouping_mat is None:
                    raise RuntimeError('Group information matrix is missing.')
                normalizer = self.sum(axis=self.Axis.HAPLOTYPE) * grouping_mat
                for hid in xrange(self.num_haplotypes):
                    self.data[hid].eliminate_zeros()  # Trying to avoid numerical problem (inf or nan)
                    self.data[hid] = np.divide(self.data[hid], normalizer)
            elif axis == self.Axis.HAPLOGROUP:  # haplotype-wise & group-wise normalization on each read
                if grouping_mat is None:
                    raise RuntimeError('Group information matrix is missing.')
                for hid in xrange(self.num_haplotypes):  # normalizer is different hap-by-hap
                    normalizer = self.data[hid] * grouping_mat  # Sparse matrix of |reads| x |loci|
                    self.data[hid].eliminate_zeros()  # Trying to avoid numerical problem (inf or nan)
                    self.data[hid] = np.divide(self.data[hid], normalizer)
            else:
                raise RuntimeError('The axis should be 0, 1, 2, or 3.')
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def pull_alignments_from(self, reads_to_use, shallow=False):
        """
        Pull out alignments of certain reads

        :param reads_to_use: numpy array of dtype=bool specifying which reads to use
        :param shallow: whether to copy sparse 3D matrix only or not
        :return: a new AlignmentPropertyMatrix object that particular reads are
        """
        new_alnmat = self.copy(shallow=shallow)
        for hid in xrange(self.num_haplotypes):
            hdata = new_alnmat.data[hid]
            hdata.data *= reads_to_use[hdata.indices]
            hdata.eliminate_zeros()
        if new_alnmat.count is not None:
            if type(new_alnmat.count) == csc_matrix:
                new_alnmat.count.data *= reads_to_use[new_alnmat.count.indices]
                new_alnmat.count.eliminate_zeros()
            elif type(new_alnmat.count) == np.ndarray:
                new_alnmat.count[np.logical_not(reads_to_use)] = 0
            else:
                raise RuntimeError('APM count should be either scipy.sparse.csc_matrix or numpy.ndarray')
        return new_alnmat

    def get_unique_reads(self, ignore_haplotype=False, shallow=False):
        """
        Pull out alignments of uniquely-aligning reads

        :param ignore_haplotype: whether to regard allelic multiread as uniquely-aligning read
        :param shallow: whether to copy sparse 3D matrix only or not
        :return: a new AlignmentPropertyMatrix object that particular reads are
        """
        if self.finalized:
            if ignore_haplotype:
                summat = self.sum(axis=self.Axis.HAPLOTYPE)
                nnz_per_read = np.diff(summat.tocsr().indptr)
                unique_reads = np.logical_and(nnz_per_read > 0, nnz_per_read < 2)
            else:  # allelic multireads should be removed
                alncnt_per_read = self.sum(axis=self.Axis.LOCUS).sum(axis=self.Axis.HAPLOTYPE)
                unique_reads = np.logical_and(alncnt_per_read > 0, alncnt_per_read < 2)
            return self.pull_alignments_from(unique_reads, shallow=shallow)
        else:
            raise RuntimeError('The matrix is not finalized.')

    def count_unique_reads(self, ignore_haplotype=False):
        if self.finalized:
            unique_reads = self.get_unique_reads(ignore_haplotype=ignore_haplotype, shallow=True)
            if ignore_haplotype:
                numaln_per_read = unique_reads.sum(axis=self.Axis.HAPLOTYPE)
                if self.count is None:
                    numaln_per_read.data = np.ones(numaln_per_read.nnz)
                else:
                    numaln_per_read.data = self.count[numaln_per_read.indices]
                return numaln_per_read.sum(axis=0).A.ravel()  # An array of size |num_loci|
            else:
                return unique_reads.sum(axis=self.Axis.READ)  # An array of size |num_haplotypes|x|num_loci|
        else:
            raise RuntimeError('The matrix is not finalized.')

    def count_alignments(self):
        if self.finalized:
            return self.sum(axis=self.Axis.READ)
        else:
            raise RuntimeError('The matrix is not finalized.')

    def report_alignment_counts(self, filename):
        alignment_counts = self.count_alignments()
        allelic_unique_counts = self.count_unique_reads(ignore_haplotype=False)
        locus_unique_counts = self.count_unique_reads(ignore_haplotype=True)
        cntdata = np.vstack((alignment_counts, allelic_unique_counts))
        cntdata = np.vstack((cntdata, locus_unique_counts))
        fhout = open(filename, 'w')
        fhout.write("locus\t" + "\t".join(['aln_%s' % h for h in self.hname]) + "\t")
        fhout.write("\t".join(['uniq_%s' % h for h in self.hname]) + "\t")
        fhout.write("locus_uniq" + "\n")
        for locus_id in xrange(self.num_loci):
            fhout.write("\t".join([self.lname[locus_id]] + map(str, cntdata[:, locus_id].ravel())) + "\n")
        fhout.close()

    def combine(self, other, shallow=False):
        if self.finalized and other.finalized:
            dmat = Sparse3DMatrix.combine(self, other)
            dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape
            if self.count is not None and other.count is not None:
                dmat.count = np.concatenate((self.count, other.count))
            if self.lengths is not None:
                dmat.lengths = copy.copy(self.lengths)
            if not shallow:
                dmat.hname = self.hname
                dmat.lname = copy.copy(self.lname)
                dmat.rname = np.concatenate((self.rname, other.rname))
                dmat.lid   = copy.copy(self.lid)
                dmat.rid   = dict(zip(dmat.rname, np.arange(dmat.num_reads)))
                dmat.__copy_group_info(self)
            return dmat
        else:
            raise RuntimeError('Both matrices must be finalized.')

    def apply_genotypes(self, gt_file):
        hid = dict(zip(self.hname, np.arange(self.num_haplotypes)))
        gid = dict(zip(self.gname, np.arange(len(self.gname))))
        gtmask = np.zeros((self.num_haplotypes, self.num_loci))
        with open(gt_file) as fh:
            if self.groups is not None:
                for curline in dropwhile(is_comment, fh):
                    item = curline.rstrip().split("\t")
                    g, gt = item[:2]
                    hid2set = np.array([hid[c] for c in gt])
                    tid2set = np.array(self.groups[gid[g]])
                    I, J = np.meshgrid(hid2set, tid2set, indexing='ij')
                    gtmask[I, J] = 1.0
            else:
                for curline in dropwhile(is_comment, fh):
                    item = curline.rstrip().split("\t")
                    g, gt = item[:2]
                    hid2set = np.array([hid[c] for c in gt])
                    I, J = np.meshgrid(hid2set, gid[g], indexing='ij')
                    gtmask[I, J] = 1.0
        self.multiply(gtmask, axis=2)
        for h in xrange(self.num_haplotypes):
            self.data[h].eliminate_zeros()

    def save(self, h5file, title=None, index_dtype='uint32', data_dtype=float, incidence_only=True, complib='zlib', shallow=False):
        Sparse3DMatrix.save(self, h5file=h5file, title=title, index_dtype=index_dtype, data_dtype=data_dtype, incidence_only=incidence_only, complib=complib)
        h5fh = tables.open_file(h5file, 'a')
        fil  = tables.Filters(complevel=1, complib=complib)
        if self.lengths is not None:
            h5fh.create_carray(h5fh.root, 'lengths', obj=self.lengths, title='Transcript Lengths', filters=fil)
        if self.count is not None:
            if len(self.count.shape) == 1:  # count is a vector
                h5fh.create_carray(h5fh.root, 'count', obj=self.count, title='Equivalence Class Counts', filters=fil)
            elif len(self.count.shape) == 2:  # count is 2-dim matrix
                cgroup = h5fh.create_group(h5fh.root, 'count', 'Sparse matrix components for N matrix')
                h5fh.create_carray(cgroup, 'indptr', obj=self.count.indptr.astype(index_dtype), filters=fil)
                h5fh.create_carray(cgroup, 'indices', obj=self.count.indices.astype(index_dtype), filters=fil)
                h5fh.create_carray(cgroup, 'data', obj=self.count.data.astype(index_dtype), filters=fil)
        if not shallow:
            h5fh.set_node_attr(h5fh.root, 'hname', self.hname)
            h5fh.create_carray(h5fh.root, 'lname', obj=self.lname, title='Locus Names', filters=fil)
            if self.rname is not None:
                h5fh.create_carray(h5fh.root, 'rname', obj=self.rname, title='Read Names', filters=fil)
            if self.sname is not None:
                h5fh.create_carray(h5fh.root, 'sname', obj=self.sname, title='Sample Names', filters=fil)
        h5fh.flush()
        h5fh.close()

    def get_read_data(self, rid):
        return self.get_cross_section(index=rid, axis=self.Axis.READ)

    def print_read(self, rid):
        """
        Prints nonzero rows of the read wanted

        """
        if self.rname is not None:
            print(self.rname[rid])
            print('--')
        r = self.get_read_data(rid)
        aligned_loci = np.unique(r.nonzero()[1])
        for locus in aligned_loci:
            nzvec = r[:, locus].todense().transpose()[0].A.flatten()
            if self.lname is not None:
                print(self.lname[locus])
            else:
                print(locus)
            print(nzvec)

    #
    # For future use
    #
    def get_reads_aligned_to_locus(self, lid, hid=None):
        ridset = set()
        if hid is None:
            for hid in xrange(self.num_haplotypes):
                curset = set(np.nonzero(self.data[hid][:, lid])[0])
                ridset = ridset.union(curset)
            return sorted(list(ridset))
        else:
            return sorted(np.nonzero(self.data[hid][:, lid])[0])


if __name__ == "__main__":
    pass  # TODO: Put a simple usage example here
