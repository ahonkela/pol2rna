import numpy as np
import scipy.io
import os

POL2DATAPATH = os.path.expanduser('~/work/projects/pol2rnaseq/')

def load_data():
    bindata = scipy.io.loadmat(POL2DATAPATH + 'all_gene_pol2bins_2013_01_02.mat')
    bininfo = scipy.io.loadmat(POL2DATAPATH + 'bininfo_dec2012_corrected.mat')
    genes = ['ENSG%011d' % x for x in bininfo['bininfo'][:,4]]
    pol2bins = bindata['pol2bins']
    d = dict()
    for k in range(len(genes)):
        if pol2bins[k,0].shape[1] > 0:
            d[genes[k]] = np.hstack(pol2bins[k,:])
    return d
