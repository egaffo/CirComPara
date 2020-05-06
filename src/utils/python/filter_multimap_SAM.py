#!/usr/bin/env python

'''
Filter SAM alignments by the number of multple hits of each read
@author: Enrico Gaffo
@email: enrico.gaffo@unipd.it
@license: any use of this script must be approved by Enrico Gaffo
'''

from __future__ import print_function
import argparse, pysam, gzip, sys

def eprint(*args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)

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
        try:
            if read.has_tag('NH'):
                if read.get_tag('NH') <= args.max_multimaps:
                    outfile.write(read)
            ## the following should work for BWA-MEM
            ## AS  Alignment score
            ## XA  Alternative hits; format: (chr,pos,CIGAR,NM;)*
            ## XS  Suboptimal alignment score
            elif read.get_tag('AS') > 0:
                ## check if it is the main alignment
                if read.get_tag('AS') > read.get_tag('XS'):
                    ## unique alignment: secondary aligmnets are suboptimal
                    ## output only the best hit
                    outfile.write(read)
                elif read.get_tag('AS') == read.get_tag('XS'):
                    ## multiple alignments: check number of hits
                    if read.is_supplementary or \
                       (not read.has_tag('XA')) or \
                       len(read.get_tag('XA').split(';')) - 1 <= args.max_multimaps:
                        outfile.write(read)
        except KeyError as ke:
            eprint(str(ke) + str(read))
