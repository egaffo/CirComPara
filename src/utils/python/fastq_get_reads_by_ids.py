#!/usr/bin/env python

import argparse
import HTSeq
import sys

def get_reads_by_ids(reads, ids):
    ''' **reads must be a HTSeq.FastqReader object, **ids is simply a set of \
        strings representing the read ids you want to get. \
        Returns a list of HTSeq.SequenceWithQualities objects.\
    '''

    out = []
    for read in reads:
      
        if read.name.split()[0].rstrip() in ids:
            out.append(read)
    
    return out


if __name__ == '__main__':

    ## Handle parameter arguments
    parser = argparse.ArgumentParser(description="Handle FASTQ files: get \
            reads in FASTQ format from a file of reads given the list of read \
            ids you want.")

    parser.add_argument('-f', '--fastq-file', action='store', dest="fastq_file",\
            required=True, type=file, help='the FASTQ file containing the reads\
            to select from')

    parser.add_argument('-l', '--ids-list-file', action='store',\
            dest="ids_list_file", required=True, type=file, help="a file \
            containing the list of read ids to select, one id per line.")

    parser.add_argument('-e','--quality-encoding', action='store',\
            dest="quality_encoding", required=False, default='phred',\
            type=str, choices=['solexa', 'solexa-old', 'phred'], help="the \
            quality encoding used in the FASTQ file")

    args = parser.parse_args()

    ## Program logic
    fastq_file = HTSeq.FastqReader(args.fastq_file, args.quality_encoding)

    ids = set([item.rstrip() for item in args.ids_list_file.readlines()])
    
    reads = get_reads_by_ids(fastq_file, ids)

    for r in reads:
        r.write_to_fastq_file(sys.stdout)
