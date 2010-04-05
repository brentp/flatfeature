r"""
flatfeature
===========
simple, stupid, flat format for genomic features.

::

    >>> flat = Flat('data/thaliana_v8.flat', 'data/thaliana_v8.fasta')
    >>> a = flat.accn('AT1G01370')
    >>> a
    (41, '1', 'AT1G01370', 143564, 145684, '+', 'CDS', [(143773, 143824), (143773, 143824)])

    >>> Flat.row_to_dict(a)
    {'ftype': 'CDS', 'accn': 'AT1G01370', 'end': 145684, 'locs': [(143773, 143824), (143773, 143824)], 'start': 143564, 'seqid': '1', 'id': 41, 'strand': '+'}

    >>> seq = flat.row_sequence('AT1G01370') 
    >>> seq == flat.row_sequence(flat[flat['accn'] == 'AT1G01370'][0])
    True

    >>> cds_seq = flat.row_cds_sequence('AT1G01370')
    >>> cds_seq == flat.row_cds_sequence(flat.accn('AT1G01370'))
    True

    >>> cds_seq[:60]
    'ATGGCGAGAACCAAGCATCGCGTTACCAGGTCACAACCTCGGAATCAAACTGATGGCGAG'
    >>> cds_seq[-60:]
    'TCAAACTGATGGCGAGAACCAAGCATCGCGTTACCAGGTCACAACCTCGGAATCAAACTG'
    >>> len(cds_seq)
    104


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

Bed
===

Bed is a subclass of Flat that provides exactly the same programmatic 
interface, but uses .bed files for storage. This is the recommended
way to use flatfeature as it is a standard format.

    >>> b = Bed('data/brachy_v1.bed.short')
    >>> bb = b.accn('Bradi1g00200')
    >>> bb
    ('Bd1', 10581, 11638, 'Bradi1g00200', '1057', '+', [(10581, 10850), (11252, 11638)], '.\t.', '.')

    >>> Bed.row_to_dict(bb)
    {'accn': 'Bradi1g00200', 'end': 11638, 'score': '1057', 'locs': [(10581, 10850), (11252, 11638)], 'start': 10581, 'rgb': '.', 'seqid': 'Bd1', 'thick': '.\t.', 'strand': '+'}

    >>> b.seqids[:4]
    ['Bd1', 'Bd5', 'scaffold_119', 'scaffold_12']

    >>> Bed.row_string(bb)
    'Bd1\t10580\t11638\tBradi1g00200\t1057\t+\t.\t.\t.\t2\t270,387\t0,671'

    >>> Bed.row_string(bb, full=False)
    'Bd1\t10580\t11638\tBradi1g00200'

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
  
def pairs_to_slice(pairs):
    """
    given a list of tuples (like a list of CDS start, stops), return
    the numpy array that will work as a slice for those tuples
    """
    return np.concatenate([np.arange(s0-1, s1) for s0, s1 in pairs])
 

class Flat(np.ndarray):
    names = ('id', 'seqid', 'accn', 'start', 'end', 'strand', 'ftype', 'locs')
    formats = ('i4', 'S32', 'S64', 'i4', 'i4', 'S1', 'S32', 'O')
    def __new__(cls, path, fasta_path=None):
        obj = np.loadtxt(path, delimiter="\t", dtype={'names': cls.names, 
                                                      'formats': cls.formats},
                                         skiprows=1, converters={7: _loc_conv})
        obj = obj.view(cls)
        obj.path = obj.filename = path
        obj.d = None
        if fasta_path is not None:
            obj.fasta = Fasta(fasta_path, flatten_inplace=True)
        return obj.view(cls)

    def __array_finalize__(self, obj):
        if obj is None: return
        self.path = self.filename = getattr(obj, 'path', getattr(obj, 'filename', None))
        self.fasta = getattr(obj, 'fasta', None)
        self.d = getattr(obj, 'd', None)

    __array_priority__ = 10

    @classmethod
    def row_to_dict(self, row):
        return dict((name, row[name]) for name in self.names)

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
        return "".join(fa[start - 1:end] for start, end in locs)

    def fill_dict(self, force=False):
        """ allow fast lookup of a row by accn name"""
        if self.d is None or force:
            self.d = dict((row['accn'], row) for row in self)


    def accn(self, accn, first_only=True):
        r = self[self['accn'] == accn]
        if first_only: 
            try:
                return r[0]
            except:
                raise KeyError("%s not found" % (accn,))
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
            fa = np.array(self.fasta[seqid]).copy()
            s = fa.shape[0]
            for row in self[self['seqid'] == seqid]:
                for start, end in row['locs']:
                    fa[start - 1: end] = mask_with
            assert s == fa.shape[0]
            if out is None:
                yield seqid, fa
            else:
                print >> out, ">%s\n%s" % (seqid, fa.tostring())

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
        raise StopIteration

    def genic_fasta(self, outfile=sys.stdout, header_key='accn'):
        if outfile: # force through iteration of generator.
            return list(self._fasta(outfile, self.row_sequence, header_key))
        else:
            return self._fasta(outfile, self.row_sequence, header_key)
    def cds_fasta(self, outfile=sys.stdout, header_key='accn'):
        if outfile: # force through iteration of generator.
            return list(self._fasta(outfile, self.row_cds_sequence, header_key))
        else:
            return self._fasta(outfile, self.row_cds_sequence, header_key)


class Bed(Flat):
    names = ('seqid', 'start', 'end', 'accn', 'score', 'strand', 'locs', 'thick', 'rgb')
    formats = ('S32', 'i4', 'i4', 'S64', 'S10', 'S1', 'O', 'S10', 'S10')
    def __new__(cls, path, fasta_path=None):
        a = []
        for line in open(path):
            if line[0] == "#": continue
            line = line.split("\t")
            L = len(line)
            if L < 12:
                line.extend(["."] * (12 - L) )

            start = int(line[1]) + 1
            end = int(line[2])
            locs = [(start, end)]
            if L == 12:
                lens = map(int, line[10].split(","))
                rel_starts = map(int, line[11].split(","))
                starts = [start + rs for rs in rel_starts]

                ends = [starts[i] + lens[i] - 1 for i in range(len(starts))]
                locs = zip(starts, ends)
            #         seqid,          end,        accn,
            a.append((line[0], start, int(line[2]), line[3], 
                      # score, strand,        # thicks 
                      line[4], line[5], locs, line[6] + "\t" + line[7],
                      line[8]))

        obj = np.array(a, dtype=zip(Bed.names, Bed.formats))
        obj = obj.view(cls)
        obj.path = obj.filename = path
        obj.d = None
        
        if fasta_path is not None:
            obj.fasta = Fasta(fasta_path, flatten_inplace=True)
        return obj.view(cls)

    @classmethod
    def row_string(cls, row, full=True):
        if not full:
            return "\t".join((row['seqid'], str(row['start'] - 1), str(row['end']), row['accn']))
        starts = [s[0] - 1 for s in row['locs']]
        ends = [s[1] for s in row['locs']]
        slens = ",".join([str(e - s) for s, e in zip(starts, ends)])
        sstarts = ",".join("%i" % (s - row['start'] + 1) for s in starts)
        return "\t".join(map(str, [row['seqid'], row['start'] - 1, row['end'], 
                                   row['accn'], row['score'], row['strand'], row['rgb'], row['thick'],
                                   len(row['locs']), slens, sstarts]))




if __name__ == "__main__":
    import doctest
    doctest.testmod()

