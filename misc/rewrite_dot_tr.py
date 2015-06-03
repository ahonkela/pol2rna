#! /usr/bin/env python

import os
import re

MYDIR=os.path.expanduser('~/local/projects/synergy_rnaseq/')
INFILE=MYDIR+'Homo_sapiens.GRCh37.68.cdna.new_ref.tr'
OUTFILE=MYDIR+'Homo_sapiens.GRCh37.68.cdna.new_ref.fixed.tr'

UNSPLICERE=re.compile(r""".*_unspliced.*""")

with open(INFILE, 'r') as fi:
    with open(OUTFILE, 'w') as fo:
        for l in fi:
            if UNSPLICERE.match(l):
                t = l.split()
                t[0] += '_unspliced'
                fo.write(' '.join(t) + '\n')
            else:
                fo.write(l)
