import os
import re
import numpy as np
import math
import cPickle

DATAPATH="/cs/taatto/group/mlb/synergy/data/2012-03_RNA"
PREMRNA_RE=re.compile(r""".*unspliced""")
READLEN=35.0
BINWIDTH=200

def read_wig(fname, matchre):
    wigs = dict()
    with open(fname) as f:
        l = f.next()
        t = l.strip().split(' ')
        curseq = t[1].split('=')[1]
        if matchre.match(curseq):
            recordme = True
            curpos = int(t[2].split('=')[1])
            curwig = np.zeros(10000000/BINWIDTH, dtype=float)
        else:
            recordme = False
        for l in f:
            if l.startswith('fixedStep'):
                t = l.strip().split(' ')
                newseq = t[1].split('=')[1]
                if matchre.match(newseq):
                    if newseq != curseq:
                        if len(wigs) % 100 == 0:
                            print len(wigs), "sequences read"
                        if newseq in wigs:
                            print "Previous record found (!)"
                        if recordme:
                            curwig.resize(curpos/BINWIDTH+1)
                            wigs[curseq] = curwig
                        curwig = np.zeros(10000000/BINWIDTH, dtype=float)
                        curseq = newseq
                    curpos = int(t[2].split('=')[1])
                    recordme = True
                else:
                    recordme = False
            else:
                if recordme:
                    curwig[curpos/BINWIDTH] += float(l)/READLEN
                    curpos += 1
    return wigs


f = os.listdir(DATAPATH)
wigs = [DATAPATH + '/' + x for x in f if x.endswith('.wig')]

w = [(x, read_wig(x, PREMRNA_RE)) for x in wigs]
with open('MCF7_wiggles.pyd', 'wb') as f:
    cPickle.dump(w, f, 2)
