import numpy as np
import scipy.io
import os

POL2DATAPATH = os.path.expanduser('~/projects/pol2rnaseq/data/')

# def summarise_halves(d):
#     return {k: (sum(v[(v.shape[0]/2):], 0)) / (sum(v, 0) + 1)
#             for k, v in d.items() if v.shape[0] > 1}

# def summarise_halves2(d):
#     return {k: np.mean(v[5:8]) - np.mean(v[2:5]) for k,v in d.items()}

def rebin(x, nbins):
    inbins = x.shape[0]
    v = np.zeros((nbins, x.shape[1]))
    bounds = np.linspace(0, inbins, nbins+1)
    for k in range(nbins):
        done = False
        #print bounds[k], '->', bounds[k+1]
        if bounds[k] < np.ceil(bounds[k]): # start with an incomplete bin
            if bounds[k+1] > np.ceil(bounds[k]):  # continue until next bin
                #print '  ', bounds[k], '->', np.ceil(bounds[k])
                v[k,] += (np.ceil(bounds[k]) - bounds[k]) * x[np.floor(bounds[k]),]
            else:   # start and finish within the same bin
                v[k,] += (bounds[k+1] - bounds[k]) * x[np.floor(bounds[k]),]
                done = True
        if not done:
            #print '  ', np.ceil(bounds[k]), '->', np.floor(bounds[k+1])
            v[k,] += np.sum(x[np.ceil(bounds[k]):(np.floor(bounds[k+1])),], 0)
            if bounds[k+1] > np.floor(bounds[k+1]):
                #print '  ', np.floor(bounds[k+1]), '->', bounds[k+1]
                v[k,] += (bounds[k+1] - np.floor(bounds[k+1])) * x[np.floor(bounds[k+1]),]
    return v

def splitratio(x, k1, k2=None):
    if k2 is not None:
        return (np.sum(x[0:k1,], 0) / (np.sum(x[0:k2,], 0) + 1e-100))
    else:
        return (np.sum(x[0:k1,], 0) / (np.sum(x, 0) + 1e-100))

def summarise_last_nth(d, n, n2=None):
    return {k: splitratio(rebin(v, n), 1, n2)
            for k, v in d.items() if v.shape[0] > 1}

def save_summaries(fits2, fname):
    with open(fname, 'w') as f:
        for k, v in sorted(fits2.items()):
            f.write('%s\t%s\n' % (k, '\t'.join(['%f' % x for x in v])))


def load_data():
    bindata = scipy.io.loadmat(POL2DATAPATH + 'all_gene_pol2bins_2014_11_19.mat')
    bininfo = scipy.io.loadmat(POL2DATAPATH + 'bininfo_nov2014_corrected.mat')
    genes = ['ENSG%011d' % x for x in bininfo['bininfo'][:,4]]
    pol2bins = bindata['pol2bins']
    d = dict()
    for k in range(len(genes)):
        if pol2bins[k,0].shape[1] > 0:
            d[genes[k]] = np.hstack(pol2bins[k,:])
    return d


def load_groseq_data():
    bindata = scipy.io.loadmat(POL2DATAPATH + 'allprof.mat')
    genes = ['ENSG%011d' % x for x in bindata['finalgenes'][:,0]]
    grobins = bindata['allprof']
    geneinfo = bindata['finalgenes']
    d = dict()
    for k in range(len(genes)):
        if geneinfo[k,3] == 1:  ## check strand
            d[genes[k]] = np.hstack(np.abs(grobins[k,:,1]))
        else:
            d[genes[k]] = np.hstack(np.abs(grobins[k,:,0]))
    return d

def analyse_pol2(pol2data):
    # pol2data = load_data()
    pol2halves = summarise_halves(pol2data)
    pol2halves2 = summarise_halves2(pol2halves)
    save_fits2(pol2halves2, 'pol2_halfdiff_2015-06-10.txt')
    return pol2halves, pol2halves2

def analyse_last5pct(dd):
    # pol2data = load_data()
    d = summarise_last_nth(dd, 20)
    save_summaries(d, 'pol2_last5pct_2015-06-10.txt')
    d = summarise_last_nth(dd, 20, 10)
    save_summaries(d, 'pol2_last5pct_b_2015-06-10.txt')
    d = summarise_last_nth(dd, 20, 2)
    save_summaries(d, 'pol2_last5pct_c_2015-06-10.txt')
