#!/usr/bin/env python
import pysam
from collections import defaultdict
from random import random
import sys
import os

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


if len(sys.argv) == 5:
    
    input = sys.argv[1]
    proportion = sys.argv[2]
    output1 = sys.argv[3]
    output2 = sys.argv[4]
    inbam = pysam.AlignmentFile(input, 'rb')
    outbam1 = pysam.AlignmentFile(output1, 'wb', template=inbam)
    outbam2 = pysam.AlignmentFile(output2, 'wb', template=inbam)

    for read1, read2 in read_pair_generator(inbam):
        p = random()
        if p <= float(proportion):
            outbam1.write(read1)
            outbam1.write(read2)
        else:
            outbam2.write(read1)
            outbam2.write(read2)

    outbam1.close()
    outbam2.close()
    inbam.close()

    pysam.sort('-o','tmp.bam',output1)
    os.rename('tmp.bam',output1)
    pysam.index(output1) 
    pysam.sort('-o','tmp.bam',output2)
    os.rename('tmp.bam',output2)
    pysam.index(output2) 
else:
    sys.exit("usage: %s <bam> <proportion> <outbam1> <outbam2> \nexample: %s /home/abc/input.bam 0.5 /home/abc/output1.bam /home/abc/output2.bam" % (sys.argv[0], sys.argv[0]))
