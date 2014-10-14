import os, sys
import numpy as np
import pysam

DATAPATH="/cs/taatto/group/mlb/synergy/data/2012-03_RNA"

TRANSCRIPTS = ["ENST00000561552", "ENST00000361711", "ENST00000343939", "ENST00000393306", "ENST00000565123", "ENST00000563543", "ENST00000567707", "ENSG00000140961_unspliced"]

UNSPLICED = ["ENSG00000140961_unspliced"]

def find_files():
    f = os.listdir(DATAPATH)
    return sorted([DATAPATH + '/' + x for x in f if x.startswith('MCF7_t')])

def gettids(f, seqs):
    return [f.gettid(x) for x in seqs]

def count_alignments(f, mytids, logf=sys.stdout):
    curqname = "invalid"
    hasother = np.zeros(len(mytids), np.bool)
    hascorrect = np.zeros(len(mytids), np.bool)
    count = np.zeros(len(mytids), np.int64)
    lines_done = 0
    for a in f:
        lines_done += 1
        if lines_done % 100000 == 0:
            logf.write(str(lines_done) +  " lines done, count: " + str(count) + "\n")
            logf.flush()
        # New alignment for the same read
        if a.qname == curqname:
            match = np.array([a.tid in x for x in mytids], np.bool)
            hascorrect = hascorrect | match
            hasother = hasother | (~match)
        # new read
        else:
            count += (hascorrect & ~hasother)
            match = np.array([a.tid in x for x in mytids], np.bool)
            hascorrect = match
            hasother = ~match
            curqname = a.qname
    return count


if __name__ == "__main__":
    myind = int(sys.argv[1])
    if len(sys.argv) > 2:
        logf = open(sys.argv[2], 'w')
    else:
        logf = sys.stdout
    fn = find_files()
    logf.write("Processing file " + fn[myind] + "\n")
    f = pysam.Samfile(fn[myind], "rb" )
    alltids = set(gettids(f, TRANSCRIPTS))
    premrnatid = set(gettids(f, UNSPLICED))
    mycount = count_alignments(f, [alltids, premrnatid], logf)
    f.close()
    f = open(DATAPATH + '/OSGIN1_counts_%d.txt' % myind, 'w')
    f.write(str(mycount) + '\n')
    f.close()
    if len(sys.argv) > 2:
        logf.close()
