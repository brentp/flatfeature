"""
>>> f = Fat("data/thaliana_v9_genes.gff.short")
>>> sorted(f.iterkeys())[:2]
['AT1G01010', 'AT1G01020']

>>> f.seqids
['1']

>>> g = f['AT1G01020']
>>> g.splicings.keys()
['.1', '.2', '.1-Protein', '.2-Protein', 'AT1G01020']

>>> g.types
['CDS', 'exon', 'five_prime_UTR', 'gene', 'mRNA', 'protein', 'three_prime_UTR']

>>> g.splicings['.1']['CDS']
[(6915L, 7069L), (7157L, 7232L), (7384L, 7450L), (7564L, 7649L), (7762L, 7835L), (7942L, 7987L), (8236L, 8325L), (8417L, 8464L), (8571L, 8666L)]

>>> g.CDS # merged of all splicings.
[[6915L, 7069L], [7157L, 7232L], [7315L, 7450L], [7564L, 7649L], [7762L, 7835L], [7942L, 7987L], [8236L, 8325L], [8417L, 8464L], [8571L, 8666L]]

>>> g.start, g.end, g.accn
(5928L, 8737L, 'AT1G01020')

>>> g = f['AT1G01050']
>>> g
Accn('AT1G01050', -, 31170:33153)

>>> list(f.get_features_in_region('1', 6915, 14127))
[Accn('AT1G01020', -, 5928:8737), Accn('AT1G01030', -, 11649:13714)]


>>> list(f.upstream(g, 1000, noncoding=False)) # features
[Accn('AT1G01060', -, 33379:37840)]

>>> list(f.downstream(g, 1000, noncoding=False)) # features
[Accn('AT1G01040', +, 23146:31227)]


>>> f.downstream(g, 12000, noncoding=True) # doctest:+ELLIPSIS
[(19169L, 23518L), (24452L, 24541L), (24656L, 24751L), (24963L, 25040L), (25436L, 25523L), (25744L, 25824L), (25998L, 26080L), (26204L, 26291L), (26453L, 26542L), (26777L, 26861L), (27013L, 27098L), (27282L, 27371L), (27534L, 27617L), (27714L, 27802L), (28432L, 28707L), (28806L, 28889L), (29081L, 29159L), (30066L, 30146L), (30312L, 30409L), (30817L, 30901L), (31080L, 31169L)]

>>> g = f['AT1G01050']
>>> f.introns(g)
[(31170L, 31381L), (31425L, 31520L), (31603L, 31692L), (31814L, 31932L), (31999L, 32087L), (32196L, 32281L), (32348L, 32430L), (32460L, 32546L), (32671L, 33153L)]


"""
import gt
gt.warning_disable()
import operator
import numpy as np
import sys
import itertools
from_iterable = itertools.chain.from_iterable

def flatten(nli):
    return list(from_iterable(nli))

class Accn(object):

    def __repr__(self):
        return "Accn('%s', %s, %i:%i)" % (self.accn, self.strand, 
                                          self.start, self.end)
    
    def _get_merged_type(self, feature_type):
        """
        this merges overlapping CDSs or mRNA's or
        a set of non overlapping, by making larger
        cdss from the overlaps
        """
        cdss = []
        for splice, ftypes in self.splicings.iteritems():
            if not feature_type in ftypes: continue
            cds = ftypes[feature_type]
            cdss.extend(cds)
        if cdss == []: return []
        return Accn._merge(cdss)

    def __getitem__(self, key):
        return self.splicings[key]


    @classmethod
    def _merge(self, flist):
        if not flist: return []
        flist = sorted(list(t) for t in flist)
        newl = [flist[0]]
        for i in range(1, len(flist)):
            # if the previous end is < next start, append
            if newl[-1][1] < flist[i][0]:
                newl.append(flist[i])
            # if the previous end is >= next start, extend previous end.
            else:
                newl[-1][1] = flist[i][1]
        return newl            


def make_accn(accn, type_loc_dict):
    _self = {'seqid': type_loc_dict.pop('seqid'),
         'strand': type_loc_dict.pop('strand'),
            "__slots__": (),
            "accn": accn,
           'splicings' : {}}
    ftypes = set([])
    
    start = int(1e30)
    end = 0
    loclist = []
    for splice, d in type_loc_dict.iteritems():
        #print d
        #print [ss for ss in [loclist for loclist in d.values()]]
        ll = []
        map(ll.extend, d.values())
        loclist.extend(ll)

        assert splice.startswith(accn), (splice, accn)
        if splice != accn:
            splice = splice.replace(accn, "")
        _self['splicings'][splice] = d
        ftypes.update(d.keys())
    _self['types'] = sorted(ftypes)
    _self['start'] = min(l[0] for l in loclist)
    _self['end'] = max(l[1] for l in loclist)

    # this creates .CDS, .gene properties. could just
    # use getattr, but then dont have nice completion for 
    # ipython.
    def gen_type_fun(ftype):
        # this becomes the method.
        def inner(self):
            return self._get_merged_type(ftype)
        return inner

    for ftype in ftypes:
        _self[ftype] = property(gen_type_fun(ftype))

    return type('Accn', (Accn,), _self)()


class Fat(object):
    seqids = None
    posns = {}

    def __init__(self, gff_filename):
        self.parse_file(gff_filename)
        self._set_start_stops()

    def __getitem__(self, accn):
        return self.accns[accn]


    def _set_start_stops(self):
        """go through and save an ordered list of positions 
        by seqid in order so we can do a binary search 
        by start, end."""
        posns = {}
        slen = 0
        for name, accn in self.accns.iteritems():
            if not accn.seqid in posns: 
                posns[accn.seqid] = {'start': [], 'end':[]}
            if len(name) > slen: slen = len(name)
            posns[accn.seqid]['start'].append((accn.start, name))
            posns[accn.seqid]['end'].append((accn.end, name))

        for seqid in posns:
            for k in ('start', 'end'):
                locs = posns[seqid][k]
                locs = np.array(locs, dtype=\
                    np.dtype([(k, np.uint32), ('accn', 'S%i' % slen)]))

                locs.sort(order=(k,))
                #if k == 'end': locs = locs[::-1]
                posns[seqid][k] = locs
        self.posns = posns
    
    def get_features_in_region(self, seqid, start, end):
        posns = self.posns[seqid]
        ix0 = posns['end']['end'].searchsorted(start, "right")
        ix1 = posns['start']['start'].searchsorted(end, "left")

        # have to grab from both to ensure we get everything.
        accns_0 = posns['start'][ix0:ix1]['accn']
        accns_1 = posns['end'][ix0:ix1]['accn']
        accns = np.unique(np.concatenate((accns_0, accns_1)))
        for accn in (self.accns[str(a)] for a in accns):
            if accn.start <= end and accn.end >= start:
                yield accn



    def upstream(self, accn, bp, noncoding=True, exclude=("CDS",)):
        """
        get the 'stuff' within bp upstream of `accn`
        if `noncoding` is true, 'stuff' is a list of
        intervals excluding any genes (or introns).
        if `noncoding` is False, 'stuff' is a list of
        Accns within that distance of bp.
        """
        if accn.strand == "+":
            end = accn.start - 1 
            start = max(end - bp, 1)
        else:    
            assert accn.strand == "-"
            start = accn.end + 1
            end = start + bp
        return self._xstream(accn, start, end, noncoding, exclude)

    def introns(self, accn):
        return self._xstream(accn, accn.start, accn.end , True, ("CDS",))

    def _xstream(self, accn, start, end, noncoding, exclude):

        feats = list(self.get_features_in_region(accn.seqid, start, end))
        if noncoding:
            all = []
            for e in exclude:
                all.extend(flatten([getattr(f, e, []) for f in feats]))
            all = Accn._merge(all)

            posns = [start]
            for s0, s1 in all:
                if s1 >= end: 
                    posns.append(s0 - 1)
                    break
                if s0 < posns[0]:
                    continue

                posns.append(s0 - 1)
                posns.append(s1 + 1)
            if len(posns) == 1: posns.append(end)

            if posns[-1] < end: 
                posns.append( end)
            else:
                posns[-1] = end
            return zip(posns[0::2], posns[1::2])
        else:
            return feats

    def downstream(self, accn, bp, noncoding=True, exclude=("CDS",)):
        if accn.strand == "-":
            end = accn.start - 1 
            start = max(end - bp, 1)
        else:    
            assert accn.strand == "+"
            start = accn.end + 1
            end = start + bp
        return self._xstream(accn, start, end, noncoding, exclude)

    def iteritems(self):
        return self.accns.iteritems()

    def iterkeys(self):
        return self.accns.iterkeys()

    def itervalues(self):
        return self.accns.itervalues()

    def parse_file(self, gff_filename):
        fi = gt.FeatureIndexMemory()
        fi.add_gff3file(gff_filename)
        self.seqids = sorted(fi.get_seqids())
        accns = {}
        for seqid in self.seqids:
            for i, parent in enumerate(sorted(fi.get_features_for_seqid(seqid), 
                                       key=operator.attrgetter('start'))):

                d = {'strand': parent.strand, 'seqid': seqid}
                ids = []
                for f in parent:
                    if "ID" in f.attribs: ids.append(f.attribs["ID"])
                    if not ids[-1] in d: d[ids[-1]] = {}
                    if not f.type in d[ids[-1]]: d[ids[-1]][f.type] = []
                    d[ids[-1]][f.type].append((f.start, f.end))

                # sometimes, the stuff is sorted by start, which puts the weird
                # names first.
                ids.sort()
                accns[ids[0]] = make_accn(ids[0], d)

        self.accns = accns


if __name__ == "__main__":

    import doctest
    doctest.testmod()
    import cPickle
    #f = Fat("data/thaliana_v9_genes.gff")
    #cPickle.dump(f, open('t.pkl', 'wb'), -1)
