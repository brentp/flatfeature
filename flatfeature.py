r"""
flatfeature
===========
simple, stupid, flat format for genomic features.

::

    >>> flat = Flat('data/thaliana_v8.flat', 'data/thaliana_v8.fasta')
    >>> flat.accn('AT1G01370')
    (41, '1', 'AT1G01370', 143564, 145684, '+', 'CDS', '143773,143824,143773,143824')

    >>> seq = flat.row_sequence('AT1G01370') 
    >>> seq == flat.row_sequence(flat[flat['accn'] == 'AT1G01370'][0])
    True

    >>> cds_seq = flat.row_cds_sequence('AT1G01370')
    >>> cds_seq == flat.row_cds_sequence(flat.accn('AT1G01370'))
    True

    >>> cds_seq[:10]
    'ATGGCGAGAA'

    >>> flat.row_locs('AT1G01370')
    [(143773, 143824), (143773, 143824)]

    >>> flat.accn('AT1G01370')['locs']
    '143773,143824,143773,143824'

    >>> flat.accn('AT1G01370')
    (41, '1', 'AT1G01370', 143564, 145684, '+', 'CDS', '143773,143824,143773,143824')

    >>> list(flat[:5].genic_fasta(outfile=None))[4].split("\n")[0]
    '>AT1G01046'

    >>> list(flat[:5].genic_fasta(outfile=None, header_key='id'))[3].split("\n")[0]
    '>4'

and that id corresponds to the row number (+ 1) in the orignal array (and
flat file)

    >>> flat[4 - 1]['accn']
    'AT1G01040'

    >>> list(flat[:5].cds_fasta(outfile=None))[0].split("\n")[0]
    '>AT1G01010'

    >>> flat.row_introns('AT1G01010')
    [(3914, 3995), (4277, 4485), (4606, 4705), (5096, 5173), (5327, 5438)]

    >>> Flat.sequence_for_locs([1, 10], flat.fasta['1'])
    'CCCTAAACCC'

    >>> flat.get_features_in_region('1', 5000, 7000)['accn']
    Flat(['AT1G01010', 'AT1G01020'], 
          dtype='|S64')

    >>> flat.get_features_in_region('1', 4000, 4000)['accn'][0]
    'AT1G01010'

"""

import numpy as np
from pyfasta import Fasta
import functools
import sys
import operator

def checkrowtype(fn):
    """ decorator:
    if the thing passed is a string get the row for that accn name and use it
    to pass to the sequence functions """
    @functools.wraps(fn)
    def wrapper(cls, row):
        if isinstance(row, str):
            row = cls[cls['accn'] == row][0]
        return fn(cls, row)
    return wrapper

class Flat(np.ndarray):
    names = ('id', 'seqid', 'accn', 'start', 'end', 'strand', 'ftype', 'locs')
    formats = ('i4', 'S12', 'S64', 'i4', 'i4', 'S1', 'S32', 'S1024')
    def __new__(cls, path, fasta_path):
        obj = np.loadtxt(path, delimiter="\t", dtype={'names': cls.names, 
                                      'formats': cls.formats}, skiprows=1)
        obj = obj.view(cls)
        obj.path = path
        obj.fasta = Fasta(fasta_path, flatten_inplace=True)
        return obj.view(cls)

    def __array_finalize__(self, obj):
        if obj is None: return
        self.path = getattr(obj, 'path', None)
        self.fasta = getattr(obj, 'fasta', None)

    __array_priority__ = 10

    @checkrowtype
    def row_sequence(self, row):
        """
        given a row from this numpy array, or an accn name.
        get the sequence given by the start and end.
        """
        assert row.shape == ()
        seqid = row['seqid']
        f = self.fasta[seqid]
        return f[row['start'] - 1:row['end']]

    @checkrowtype
    def row_cds_sequence(self, row):
        """
        given a row from this numpy array, or an accn name.
        get the sequence given by the cds locs in the final
        column.
        """
        assert row.shape == ()
        seqid = row['seqid']
        fa = self.fasta[seqid]
        locs = map(int, row['locs'].split(","))
        return Flat.sequence_for_locs(locs, fa)

    @classmethod
    def sequence_for_locs(cls, locs, fa):
        seq = []
        for i in range(0, len(locs), 2):
            seq.append(fa[locs[i] - 1: locs[i + 1]])
        return "".join(seq)

    def accn(self, accn, first_only=True):
        r = self[self['accn'] == accn]
        if first_only: return r[0]
        return r

    def get_features_in_region(self, seqid, start, end):
        assert start <= end
        return self[
            (self['seqid'] == seqid) & 
            (self['start'] <= end) &
            (self['end'] >= start)
        ]

    def mask_cds(self, mask_with='N', out_fasta=None, feat_list=None):
        """
        mask the cds sequence for each row associated with this object.
        if out_fasta is given, the output will be written as a fasta
        to that file. other wise, tuples of (seqid, sequence)
        will be generated.
        if feat_list is specified, only features in that list will be
        masked. TODO
        """
        return self._mask(True, mask_with, out_fasta)

    def mask_genic(self, mask_with='N', out_fasta=None):
        """
        mask the gene sequence associated with this object.
        if out_fasta is given, the output will be written as a fasta
        to that file. other wise, tuples of (seqid, sequence)
        will be generated.
        if feat_list is specified, only features in that list will be
        masked. TODO
        """
        return self._mask(False, mask_with, out_fasta)


    def _mask(self, cds, mask_with, out):
        """ 
        yields tuples of seqid, masked_sequence
        for each chromsome in fasta
        where masked_sequence has all cds sequence
        masked.
        """
        if out is not None and isinstance(out, basestring):
            out = open(out, 'wb')
        for seqid in sorted(self.fasta.keys()):
            fa = np.array(self.fasta[seqid].copy())
            s = fa.shape[0]
            for row in self[self['seqid'] == seqid]:
                locs = self.row_locs(row) if cds else [(row['start'], row['end'])]
                for start, end in locs:
                    fa[start - 1: end] = mask_with
            assert s == fa.shape[0]
            if out is None:
                yield seqid, fa
            else:
                print >> out, ">%s\n%s" % (seqid, fa)

    def row_introns(self, row):
        """
        grap the introns for this feature
        """
        locs = reduce(operator.add, self.row_locs(row))[1:-1]

        its = zip([x + 1 for x in locs[0::2]], 
                  [x - 1 for x in locs[1::2]])
        return its

    @checkrowtype
    def row_locs(self, row):
        try:
            locs = map(int, row['locs'].split(","))
        except:
            print row
            raise
        return zip(locs[::2], locs[1::2])

    row_exons = row_locs

    def _fasta(self, outfile, seq_fn, header_key):
        if isinstance(outfile, basestring):
            outfile = open(outfile, 'w')
        for row in self:
            header = row[header_key]
            if outfile is None:
                yield ">%s\n%s" % (header, seq_fn(row))
            else:
                print >>outfile, ">%s" % header
                print >>outfile, seq_fn(row)

    def genic_fasta(self, outfile=sys.stdout, header_key='accn'):
        return self._fasta(outfile, self.row_sequence, header_key)
    def cds_fasta(self, outfile=sys.stdout, header_key='accn'):
        return self._fasta(outfile, self.row_cds_sequence, header_key)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
