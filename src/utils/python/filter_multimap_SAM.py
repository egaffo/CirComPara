#!/usr/bin/env python

'''
Filter SAM alignments by the number of multple hits of each read
@author: Enrico Gaffo
@email: enrico.gaffo@unipd.it
@license: any use of this script must be approved by Enrico Gaffo
'''

import argparse, pysam, gzip

if __name__ == '__main__':
    # just handle input parameters
    parser = argparse.ArgumentParser(description = \
                                        '''Filter SAM alignments by the number of multple '''\
                                        '''hits of each read. Output to stdout.''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('samfile', type = str, default = '-',
    	 	    help = 'Alignment file in SAM format. Can read Gzipped compressed SAM files.')
    parser.add_argument('-m', '--max_multimaps', default = float('inf'), type = int,
            help = 'Filter out read with more than MAX_MULTIMAPS multiple alignments.')
    args = parser.parse_args()
    

    infile = args.samfile
    # if compressed file use appropriate function
    if infile.endswith('.gz'):
        infile = gzip.open(infile, 'r')
    
    samfile = pysam.AlignmentFile(infile)
    
    outfile = pysam.AlignmentFile("-", "w", template = samfile)

    for read in samfile:
        if read.has_tag('NH'):
            if read.get_tag('NH') <= args.max_multimaps:
                outfile.write(read)
        # the following should work for BWA-MEM
        elif read.has_tag('XA'):
            if len(read.get_tag('XA').split(';')) - 1 <= args.max_multimaps:
                outfile.write(read)
    
