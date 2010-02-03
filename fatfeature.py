"""
>>> f = Fat('some.gff')
>>> f.keys()
['At1g10010', ... 'At1g58960']
>>> f.seqids
['1', '2', '3', '4', '5']

>>> g = f['At2g26540']
>>> g.splicings
['.1', '.2']

>>> g['.1.'].types
['CDS', 'mRNA', 'five_prime_utr', 'three_prime_utr', 'gene']

>>> g['.1'].CDS
[(100, 200), (300, 400)]

>>> g.CDS
{'.1':[(100, 200), (300, 400)], '.2': [(98, 198), (298, 398)]}

>>> Fat.merge(g.CDS)
[(98, 200), (298, 400)]

>>> g.upstream(1000, exclude='all') # only non-feature (non coding)
[(55, 99), (0, 22)] 

>>> g.upstream(1000, features=True)
['At1g...', ...]

"""
import gt
gt.warning_disable()
import operator
import pprint

class Accn(object):
    def __init__(self, type_loc_dict):
        pass        

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
                    #print subf.attribs, subf.type, subf.get_type()
                accns[ids[0]] = d
                    
                if i == 1: 
                    pprint.pprint( accns)
                    1/0
                    break
        self.accns = accns

Fat("data/thaliana_v9_genes.gff")
