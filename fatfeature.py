"""
>>> f = Fat('some.gff')
>>> f.keys()
['At1g10010', ... 'At1g58960']
>>> f.seqids
['1', '2', '3', '4', '5']

>>> g = f['At2g26540']
>>> g.splicings
['.1', '.2']

>>> g['.1.].types
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

class GFFLine(object):
    __slots__ = ('seqid', 'com', 'type', 'start', 'stop', 'end', 'strand', 'other',
                    'attrs', 'attribs', 'sattrs', 'orig')
    def __init__(self, sline):
        line = sline.rstrip().split("\t")
        self.seqid = line[0]
        self.com  = line[1]
        self.type = line[2]
        self.start = int(line[3])
        self.stop = self.end = int(line[4])
        self.orig = line[5]
        self.strand = line[6] in ('-', '-1') and -1 or 1
        self.other = line[7]
        self.sattrs = line[8]
        self.attrs = self.attribs = self._parse_attrs(line[8])

    def _parse_attrs(self, sattrs):
        attrs = {}
        if "=" in sattrs:
            for pair in sattrs.split(';'):
                if not pair: continue
                pair = pair.split('=')
                attrs[pair[0]] = pair[1]
        if attrs == {}: attrs["ID"] = sattrs.rstrip("\r\n; ")
        return attrs

    @classmethod
    def _get_non_comment_line(cls, fh):
        while True:
            line = fh.readline()
            if line and line[0] == '#': continue
            return line


    def __repr__(self):
        return "GL(%s %s %s)" % (self.seqid, self.type, self.attrs)

    
def _get_gff_block(fh):
    """
    so this just continues to append lines to lines[]
    until Parent attribute of a line does not match
    any of the values in parent_ids
    any time a line is found with a new ID attribute
    (whose Parent attr matches the current parent_ids list),
    that line's own ID is added to the parent_ids list"""

    lines = [GFFLine._get_non_comment_line(fh)]
    if _get_gff_block.parent_ids is None:
        _get_gff_block.parent_ids = [GFFLine(lines[0]).attribs["ID"]]

    parent_ids = _get_gff_block.parent_ids

    while True:
        lines.append(GFFLine._get_non_comment_line(fh))
        new_parent = False
        for parent_id in parent_ids:
            if 'Parent=' + parent_id in lines[-1]:
                if not lines[-1]: 1/0; return lines # end of file
                lines.append(GFFLine._get_non_comment_line(fh))
                if 'ID=' in lines[-1]:
                    new_parent=True
                break
        if new_parent:
            parent_ids.append(GFFLine(lines[-1]).attribs["ID"])
        else:
            break
    if len(parent_ids) == 1:
        parent_ids = [GFFLine(lines[-1]).attribs["ID"]]
        print parent_ids 
    else:
        parent_ids = [parent_ids[-1]]
    
    _get_gff_block.parent_ids = parent_ids
    return [GFFLine(l) for l in lines]

_get_gff_block.parent_ids = None
_get_gff_block.lines = None



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

    def parse_file(self, gff_filename, skip=('chromosome', )):
        fh = open(gff_filename)
        while True:
            block = _get_gff_block(fh)
            if block is None: return 
            if block[0].type in skip: continue
            try:
                print "ID:", block[0].attribs["ID"]
            except:
                print block[0]
                raise
            


Fat("data/thaliana_v9_genes.gff")
