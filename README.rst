flatfeature
===========
simple, stupid, flat format for genomic features.
all information for a given feature is saved on a single line.
the format is described by the columns::

    id  chr accn    start   stop    strand  ftype   locs

where ftype is usually CDS if available, otherwise it's the
highest level feature type. e.g. 'miRNA' or 'pseudogene'...
the locs is a string containing the start,stops.

this module requires pyfasta (which is available via easy_install) and numpy >= 1.4.0
NOTE: this module now also supports the more common Bed format! The usage
is the same as with Flat except the `Bed`_ constructor is used on a bed file.

::

    >>> from flatfeature import Flat
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

    >>> list(flat[:5].genic_fasta(outfile=None))[4].split("\n")[0]
    '>AT1G01046'

    >>> list(flat[:5].genic_fasta(outfile=None, header_key='id'))[3].split("\n")[0]
    '>4'

and that id corresponds to the row number (+ 1) in the orignal array (and
flat file) ::

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


example useage: see how mean features on chromosome 4 are on the '-' strand ::

    >>> flat[(flat['seqid'] == '4') & (flat['strand'] == '-')].shape
    (2502,)


Bed
===

Bed is a subclass of Flat that provides exactly the same programmatic
interface, but uses .bed files for storage. This is the recommended
way to use flatfeature as it is a standard format::

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
