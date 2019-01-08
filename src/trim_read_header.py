#!/usr/bin/env python

import os, sys, argparse, gzip
from collections import namedtuple

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

if __name__ == '__main__':

    usage = """
            Reduce read name to ID only, discard description if present.
            Optionally, append suffix to read ID.
            """

    parser = argparse.ArgumentParser(description = usage, 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', default = '-', help = 'A FASTQ file or stream (default -)')
    parser.add_argument('-s', '--suffix', type = str, dest = 'suffix', 
                        default = '',
                        help = 'Suffix to append to read ID. F.i. "\\1" or "\\2" '\
                               'to discriminate first/second read mates')
    args = parser.parse_args()

    ## handle compressed files 
    if args.input.endswith('.gz'):
        f = gzip.open(args.input, 'rb')
    elif args.input == '-':
        f = sys.stdin
    else:
        f = open(args.input, 'r')
    
    for fq in fq_blocks(f):
        
        print fq.name.split(" ")[0] + args.suffix
        print fq.seq
        print fq.plus
        print fq.qual
 
