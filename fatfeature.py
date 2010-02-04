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


#>>> f.upstream(g, 1000, exclude='all') # only non-feature (non coding)
[(55, 99), (0, 22)] 

#>>> f.upstream(g, 1000, features=True)
['At1g...', ...]

"""
import gt
gt.warning_disable()
import operator
import pprint

class Accn(object):

    def __repr__(self):
        return "Accn('%s', %s)" % (self.accn, self.strand)
    
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

    def __init__(self, gff_filename):
        self.parse_file(gff_filename)

    def __getitem__(self, accn):
        return self.accns[accn]

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

    #f = Fat("data/thaliana_v9_genes.gff.short")
    import doctest
    doctest.testmod()
