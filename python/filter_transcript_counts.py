# python filter_transcript_counts.py < transcript_counts.txt > active_transcripts.txt

import sys

print "Gene\tTranscript\tExpression"
for l in sys.stdin:
    t = l.strip().split('\t')
    if float(t[2]) > 1.1:
        print '\t'.join(t[0:3])
