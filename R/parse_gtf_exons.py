import os.path
import numpy as np
import pandas as pd
import gzip
import re

#fname = os.path.expanduser('~/Downloads/Homo_sapiens.GRCh37.72.gtf')
fname = os.path.expanduser('~/data/genomes/Homo_sapiens.GRCh37.68.gtf.gz')

def gtf_desc_to_dict(s):
    t = s.strip().split(';')
    t = [x.strip() for x in t]
    d = {}
    for l in t[:-1]:
        i = l.index(' ')
        key = l[0:i]
        val = l[i+1:].strip('"')
        d[key] = val
    return d

def read_exons(fname):
    genes = dict()
    for l in gzip.open(fname):
        t = l.strip().split('\t')
        if t[2] == 'exon':
            d = gtf_desc_to_dict(t[8])
            d['exon_id']='EXON%s_%s' % (t[3], t[4])
            t[8] = d
            ense = '%s.%s' % (d['gene_id'], d['exon_id'])
            if d['gene_id'] in genes:
                genes[d['gene_id']].append(d['exon_id'])
            else:
                genes[d['gene_id']] = [d['exon_id']]
    return genes

def read_transcripts(fname):
    genes = dict()
    for l in gzip.open(fname):
        t = l.strip().split('\t')
        if t[2] == 'exon':
            d = gtf_desc_to_dict(t[8])
            d['exon_id']='EXON%s_%s' % (t[3], t[4])
            t[8] = d
            if d['gene_id'] in genes:
                if d['transcript_id'] in genes[d['gene_id']]:
                    genes[d['gene_id']][d['transcript_id']].append(t)
                else:
                    genes[d['gene_id']][d['transcript_id']] = [t]
            else:
                genes[d['gene_id']] = {d['transcript_id']: [t]}
    return genes

def find_skips(exons1, exons2):
    pat1 = ''.join(['T' if k in exons1 else 'F' for k in exons2])
    pat2 = ''.join(['T' if k in exons2 else 'F' for k in exons1])
    if re.match(r"""TF+T""", pat1) or re.match(r"""TF+T""", pat2):
        return True
    else:
        return False

def find_skipped_exons(mytranscripts):
    myexons = {k: [t[8]['exon_id'] for t in v] for k,v in mytranscripts.items()}
    for i, ex1 in enumerate(myexons.values()):
        for ex2 in myexons.values()[i+1:]:
            if find_skips(ex1, ex2):
                return True
    return False

if 'exons' not in locals():
    print 'Reading exons'
    exons = read_exons(fname)

if 'transcripts' not in locals():
    print 'Reading transcripts'
    transcripts = read_transcripts(fname)

print 'Finding skips'
hasskips = [find_skipped_exons(k) for k in transcripts.values()]

geneskips = pd.DataFrame(hasskips, index=transcripts.keys())
geneskips = geneskips.sort_index()
geneskips.to_csv('skipped_exons.csv')
