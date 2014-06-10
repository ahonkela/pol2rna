import os
import re
import numpy as np
import math
import cPickle
from matplotlib import pyplot as plt

SAVEPATH="/share/synergy/analyses/premrna/"
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

def read_and_save():
    f = os.listdir(DATAPATH)
    wigs = [DATAPATH + '/' + x for x in f if x.endswith('.wig')]

    w = [(x, read_wig(x, PREMRNA_RE)) for x in wigs]
    with open(SAVEPATH + 'MCF7_wiggles.pyd', 'wb') as f:
        cPickle.dump(w, f, 2)
    return w

def load_wigs():
    with open(SAVEPATH + 'MCF7_wiggles.pyd', 'r') as f:
        w = cPickle.load(f)
    w.sort(key=lambda x: x[0])
    return w

def combine_wigs(w):
    if len(w[0]) == 2:
        w = [x[1] for x in w]
    genes = set()
    val = dict()
    for k in w:
        genes = genes.union(k.keys())
    print len(genes), "genes in all"
    for g in genes:
        d = [x.get(g, np.zeros(0, dtype=float)) for x in w]
        maxlen = max([len(x) for x in d])
        v = np.zeros((maxlen, len(d)), dtype=float)
        for i in range(len(d)):
            v[0:len(d[i]), i] = d[i]
        val[g] = v
    return val

def plot_one(u, gene=None, savefile=None):
    fig, ax = plt.subplots(figsize = (6, 4))
    u = u / np.max(u)
    u = u + np.array(range(9, -1, -1))[np.newaxis, :]
    ax.plot(u)
    ax.set_yticks(np.arange(10)+0.5)
    ax.set_yticklabels(('1280', '640', '320', '160', '80', '40', '20', '10', '5', '0'))
    if gene is not None:
        ax.set_title(gene)
    fig.tight_layout()
    if savefile is not None:
        fig.savefig(savefile, format='png')
    return fig, ax

def plot_wigs(v, genes, savepath=None):
    for g in genes:
        plt.close()
        try:
            u = v[g + '_unspliced']
            if savepath is not None:
                plot_one(u, g, savepath+'/'+g+'_premrna.png')
            else:
                plot_one(u, g)
        except:
            print "Plotting gene", g, "failed (no data?)"


w = load_wigs()
v = combine_wigs(w)

PLOTPATH="/share/synergy/analyses/premrna"
with open(os.path.expanduser('~/Dropbox/projects/pol2rnaseq/interesting_genes.txt'), 'r') as f:
    genes = [x.strip() for x in f.readlines()]

plot_wigs(v, genes, PLOTPATH)
