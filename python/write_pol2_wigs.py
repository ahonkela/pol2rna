import numpy as np
import scipy.io
import os

TIMES = ['0000', '0005', '0010', '0020', '0040',
         '0080', '0160', '0320', '0640', '1280']

NORM = 36

def load_data(datapath=os.path.expanduser('~/projects/pol2rnaseq/')):
    bindata = scipy.io.loadmat(datapath + 'all_gene_pol2bins_2013_01_02.mat')
    bininfo = scipy.io.loadmat(datapath + 'bininfo_dec2012_corrected.mat')
    genes = ['ENSG%011d' % x for x in bininfo['bininfo'][:,4]]
    pol2bins = bindata['pol2bins']
    d = dict()
    for k in range(len(genes)):
        if pol2bins[k,0].shape[1] > 0:
            d[genes[k]] = np.hstack(pol2bins[k,:]/NORM)
    return d

def write_bin_wigs(w):
    genes = sorted(w.keys())
    for i in range(len(TIMES)):
        fname = 'MCF7_PolII_t%smin.wig' % TIMES[i]
        with open(fname, 'w') as f:
            f.write('track type=wiggle_0 name="mywiggle" description="mywiggle" visibility=full\n')
            for gene in genes:
                f.write('fixedStep chrom=%s start=1 step=200 span=200\n' % gene)
                f.write('\n'.join([('%f' % x).rstrip('0').rstrip('.') for x in w[gene][:,i]]) + '\n')
