import os.path
import numpy as np
import gzip

#fname = os.path.expanduser('~/Downloads/Homo_sapiens.GRCh37.72.gtf')
fname = os.path.expanduser('~/data/genomes/Homo_sapiens.GRCh37.69.gtf.gz')

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

def read_transcript_exons(fname):
    transcripts = dict()
    for l in gzip.open(fname):
        t = l.strip().split('\t')
        if t[2] == 'exon':
            d = gtf_desc_to_dict(t[8])
            t[8] = d
            enst = '%s.%s' % (d['gene_id'], d['transcript_id'])
            if enst in transcripts:
                transcripts[enst].append(t)
            else:
                transcripts[enst] = [t]
    return transcripts

def transcript_intron_lengths(exons):
    startpos = np.array([x[3] for x in exons], dtype=int)
    endpos = np.array([x[4] for x in exons], dtype=int)
    strand = exons[0][6]
    number = np.array([x[8]['exon_number'] for x in exons], dtype=int)
    assert(all(sorted(number) == number))
    if strand == '+':
        introns = startpos[1:] - endpos[0:-1] - 1
    else:
        introns = startpos[0:-1] - endpos[1:] - 1
    return introns

def transcript_intron_lengths_and_positions(exons):
    startpos = np.array([x[3] for x in exons], dtype=int)
    endpos = np.array([x[4] for x in exons], dtype=int)
    strand = exons[0][6]
    number = np.array([x[8]['exon_number'] for x in exons], dtype=int)
    assert(all(sorted(number) == number))
    if strand == '+':
        introns = startpos[1:] - endpos[0:-1] - 1
        positions = endpos[-1]
        length = endpos[-1] - startpos[0]
        exonlen5 = endpos[0] - startpos[0]
        exonlen3 = endpos[-1] - startpos[-1]
    else:
        introns = startpos[0:-1] - endpos[1:] - 1
        positions = startpos[-1]
        length = endpos[0] - startpos[-1]
        exonlen5 = endpos[-1] - startpos[-1]
        exonlen3 = endpos[0] - startpos[0]
    return (introns, positions, strand, length, exonlen5, exonlen3)

def print_intron_lengths(fname, introns):
    f = open(fname, 'w')
    for k in sorted(introns.items()):
        f.write('%s\t%s\n' % (k[0], '\t'.join([str(x) for x in k[1]])))
    f.close()

def print_intron_lengths_and_positions(fname, introns):
    f1 = open(fname + '.lengths2', 'w')
    for k in sorted(introns.items()):
        f1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (k[0], k[1][2], k[1][1], k[1][3], k[1][4], k[1][5], '\t'.join([str(x) for x in k[1][0]])))
    f1.close()

if 'transcripts' not in locals():
    transcripts = read_transcript_exons(fname)

#introns = {enst: transcript_intron_lengths(exons) for (enst, exons) in transcripts.items()}

introns_with_pos = {enst: transcript_intron_lengths_and_positions(exons) for (enst, exons) in transcripts.items()}

#print_intron_lengths('intron_lengths.txt', introns)
print_intron_lengths_and_positions('intron_lengths.txt', introns_with_pos)
