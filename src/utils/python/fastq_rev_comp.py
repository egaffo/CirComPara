#!/usr/bin/env python

import argparse, gzip, HTSeq, sys

def get_revcomp(reads):
    
    out = []
    for read in reads:
        out.append(read.get_reverse_complement())

    return out


if __name__ == '__main__':

    ## Handle parameter arguments
    parser = argparse.ArgumentParser(description = "Return to stout \
            the reverse complement of the input FASTQ reads.")

    parser.add_argument('-f', '--fastq-file', 
            action = 'store', 
            dest = "fastq_file",
            required = True, 
            type = str, 
            help = "The FASTQ file to be reverse complemented. Gzip\'d files "\
                   "and stdin are supported. Set '-' to read from stdin")

    parser.add_argument('-e','--quality-encoding', 
            action = 'store',
            dest = "quality_encoding", 
            required = False, 
            default = 'phred',
            type = str, 
            choices = ['solexa', 'solexa-old', 'phred'], 
            help = "the quality encoding used in the FASTQ file")

    args = parser.parse_args()

    ## handle compressed files 
    if args.fastq_file.endswith('.gz'):
        f = gzip.open(args.fastq_file, 'rb')
    elif args.fastq_file == '-':
        f = sys.stdin
    else:
        f = open(args.fastq_file, 'r')

    ## Program logic
    fastq_file = HTSeq.FastqReader(f, args.quality_encoding)

    reads = get_revcomp(fastq_file)

    for r in reads:
        r.write_to_fastq_file(sys.stdout)
