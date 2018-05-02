#!/usr/bin/env python

'''
Print the length of (multi)FASTA enties
@author: Enrico Gaffo
@email: enrico.gaffo@unipd.it
@license: any use of this script must be approved by Enrico Gaffo
'''

import re, argparse, gzip, sys
from Bio import SeqIO

if __name__ == '__main__':
    # just handle input parameters
    parser = argparse.ArgumentParser(description = \
                                        '''Print the (multi)FASTA length(s)''')
    parser.add_argument('infile', type = str, default = '-',
    	 	    help = 'A (multi)FASTA file')
    args = parser.parse_args()
    
    if args.infile == '-':
        f = sys.stdin
    else:
        # if compressed file use appropriate function
        if args.infile.endswith('.gz'):
            open = gzip.open
        
        f = open(args.infile, 'r')
            
    # load FASTA entries
    transcripts = SeqIO.parse(f, "fasta")
    # cicle through the entries
    for trs in transcripts:
        sys.stdout.write(trs.id + '\t' + str(len(trs.seq)) + '\n')


