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

def output_fixed(filename, suffix):

    ## handle compressed files 
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rb')
    elif filename == '-':
        f = sys.stdin
    else:
        f = open(filename, 'r')
    
    for fq in fq_blocks(f):
        
        print fq.name.split(" ")[0] + suffix
        print fq.seq
        print fq.plus
        print fq.qual

if __name__ == '__main__':

    usage = """
            Reduce read name to ID only, discard description if present.
            Optionally, append suffix to read ID. It can concatenate 
            paired-end reads and append '1' and '2' suffixes by default.
            """

    parser = argparse.ArgumentParser(description = usage, 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--mate1', default = '-', required = True,
                        help = 'A FASTQ file or stream. Gzip\'d file '\
                               'supported. Use - for stream input')
    parser.add_argument('-r', '--mate2', default = None, required = False,
                        help = 'A FASTQ file. Gzip\'d file supported')
    parser.add_argument('-s', '--suffix', type = str, dest = 'suffix', 
                        default = '',
                        help = 'Suffix to append to read ID. '\
                               'When paired-end reads, the string is '\
                               'prepended to the default "1" and "2" suffixes')
    args = parser.parse_args()

    suffix = args.suffix

    if args.mate2:
        suffix = args.suffix + '1'
    
    output_fixed(args.mate1, suffix)

    if args.mate2:
        suffix = args.suffix + '2'
        output_fixed(args.mate2, suffix)

