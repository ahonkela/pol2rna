import numpy as np
import scipy.io
import os

TIMES = ['0000', '0005', '0010', '0020', '0040',
         '0080', '0160', '0320', '0640', '1280']

CHRS = {k:str(k) for k in range(1, 23)}
CHRS[23] = 'X'
CHRS[24] = 'Y'
CHRS[25] = 'MT'
CHRS[-1] = '-1'

NORM = 200

def load_data(datapath=os.path.expanduser('~/projects/pol2rnaseq/data/')):
    bindata = scipy.io.loadmat(datapath + 'all_gene_pol2bins_2014_11_19.mat')
    bininfo = scipy.io.loadmat(datapath + 'bininfo_nov2014_corrected.mat')
    genes0 = ['ENSG%011d' % x for x in bininfo['bininfo'][:,4]]
    genes = [''] * len(genes0)
    for k in xrange(len(genes)):
        info = bininfo['bininfo'][k,:]
        genes[k] = 'ENSG%011d:chr%s:%d-%d:%c' % (info[4], CHRS[info[0]], info[1], info[2], '+' if info[5] > 0 else '-')
    lens = {k:v for k, v in zip(genes, bininfo['bininfo'][:,2] - bininfo['bininfo'][:,1])}
    pol2bins = bindata['pol2bins']
    d = dict()
    for k in range(len(genes)):
        if pol2bins[k,0].shape[1] > 0:
            d[genes[k]] = np.hstack(pol2bins[k,:]/NORM)
    return (d, lens)

def write_bin_wigs(w0):
    w, lens = w0
    genes = sorted(w.keys())
    for i in range(len(TIMES)):
        fname = 'MCF7_PolII_t%smin.wig' % TIMES[i]
        with open(fname, 'w') as f:
            f.write('track type=wiggle_0 name="mywiggle" description="mywiggle" visibility=full\n')
            for gene in genes:
                myend = lens[gene]+2
                steps = np.concatenate((np.array((1,)), np.arange(myend % 200, myend, 200)))
                f.write('variableStep chrom=%s span=200\n' % gene)
                f.write('\n'.join([('%d %f' % (y, x)).rstrip('0').rstrip('.') for x, y in zip(w[gene][::-1,i], steps)]) + '\n')
