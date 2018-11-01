#!/usr/bin/env python

import os,sys
from optparse import *
import gzip
from collections import namedtuple

usage = """

  %prog <reads.fastq.gz> > trimmed_header_reads.fastq

Reduce read name to ID only, discard description if present. Can handle gzipped files.
"""

parser = OptionParser(usage = usage)
options,args = parser.parse_args()

def fq_blocks(lines):
    """
    Iterates over FASTQ text lines and yields
    NamedTuple instances of the FASTQ with attributes 
    fq.name, fq.seq, fq.plus, fq.qual
    """
    FQ = namedtuple("fq_tuple", "name, seq, plus, qual")
    
    I = lines.__iter__()

    try:
        while I:
            name = I.next().rstrip()
            seq = I.next().rstrip()
            plus = I.next().rstrip()
            qual = I.next().rstrip()
        
            yield FQ(name, seq, plus, qual)
    except StopIteration:
        pass

## handle compressed files 
if args[0].endswith('.gz'):
    f = gzip.open(args[0], 'rb')
else:
    f = open(args[0], 'r')

for fq in fq_blocks(f):
    
    print fq.name.split(" ")[0]
    print fq.seq
    print fq.plus
    print fq.qual
 
