r"""
flatfeature
===========
simple, stupid, flat format for genomic features.

::

    >>> flat = Flat('data/thaliana_v8.flat', 'data/thaliana_v8.fasta')
    >>> flat.accn('AT1G01370')
    (41, '1', 'AT1G01370', 143564, 145684, '+', 'CDS', [(143773, 143824), (143773, 143824)])

    >>> seq = flat.row_sequence('AT1G01370') 
    >>> seq == flat.row_sequence(flat[flat['accn'] == 'AT1G01370'][0])
    True

    >>> cds_seq = flat.row_cds_sequence('AT1G01370')
    >>> cds_seq == flat.row_cds_sequence(flat.accn('AT1G01370'))
    True

    >>> cds_seq[:10]
    'ATGGCGAGAA'

    >>> flat.accn('AT1G01370')['locs']
    [(143773, 143824), (143773, 143824)]

    >>> flat.accn('AT1G01370')
    (41, '1', 'AT1G01370', 143564, 145684, '+', 'CDS', [(143773, 143824), (143773, 143824)])

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

    >>> Flat.sequence_for_locs([(1, 10)], flat.fasta['1'])
    'CCCTAAACCC'

    >>> flat.get_features_in_region('1', 5000, 7000)['accn']
    Flat(['AT1G01010', 'AT1G01020'], 
          dtype='|S64')

    >>> flat.get_features_in_region('1', 4000, 4000)['accn'][0]
    'AT1G01010'

    >>> Flat.row_string(flat[0])
    '1\t1\tAT1G01010\t3631\t5899\t+\tCDS\t3760,3913,3996,4276,4486,4605,4706,5095,5174,5326,5439,5630'

    >>> flat.seqids
    ['1', '2', '3', '4', '5']

    >>> flat.fill_dict()
    >>> flat.d["AT1G01010"]
    (1, '1', 'AT1G01010', 3631, 5899, '+', 'CDS', [(3760, 3913), (3996, 4276), (4486, 4605), (4706, 5095), (5174, 5326), (5439, 5630)])

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
        r = row
        if isinstance(row, str):
            try:
                r = cls[cls['accn'] == row][0]
            except IndexError:
                print >>sys.stderr, row
                raise
        return fn(cls, r)
    return wrapper

def _loc_conv(locstr):
    locs = map(int, locstr.split(","))
    return zip(locs[::2], locs[1::2])
    

class Flat(np.ndarray):
    names = ('id', 'seqid', 'accn', 'start', 'end', 'strand', 'ftype', 'locs')
    formats = ('i4', 'S12', 'S64', 'i4', 'i4', 'S1', 'S32', 'O')
    def __new__(cls, path, fasta_path=None):
        obj = np.loadtxt(path, delimiter="\t", dtype={'names': cls.names, 
                                                      'formats': cls.formats},
                                         skiprows=1, converters={7: _loc_conv})
        obj = obj.view(cls)
        obj.path = path
        if fasta_path is not None:
            obj.fasta = Fasta(fasta_path, flatten_inplace=True)
        return obj.view(cls)

    def __array_finalize__(self, obj):
        if obj is None: return
        self.path = getattr(obj, 'path', None)
        self.fasta = getattr(obj, 'fasta', None)

    __array_priority__ = 10

    @classmethod
    def row_to_dict(self, row):
        return dict((name, row[name]) for name in Flat.names)

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

    @property
    def seqids(self):
        return sorted(np.unique(self['seqid']).tolist())

    @classmethod
    def row_string(cls, row):
        strlocs = ",".join("%i,%i" % pair for pair in row['locs'])
        return "\t".join(map(str, [row[c] for c in cls.names[:-1]])) + "\t" + strlocs

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
        return Flat.sequence_for_locs(row['locs'], fa)

    @classmethod
    def sequence_for_locs(cls, locs, fa):
        seq = []
        for start, end in locs:
            seq.append(fa[start - 1: end])
        return "".join(seq)

    def fill_dict(self):
        """ allow fast lookup of a row by accn name"""
        self.d = dict((row['accn'], row) for row in self)


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
                for start, end in row['locs']:
                    fa[start - 1: end] = mask_with
            assert s == fa.shape[0]
            if out is None:
                yield seqid, fa
            else:
                print >> out, ">%s\n%s" % (seqid, fa)

    @checkrowtype
    def row_introns(self, row):
        """
        grap the introns for this feature
        """
        locs = reduce(operator.add, row['locs'])[1:-1]

        its = zip([x + 1 for x in locs[0::2]], 
                  [x - 1 for x in locs[1::2]])
        return its

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
