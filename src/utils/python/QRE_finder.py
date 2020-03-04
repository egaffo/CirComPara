#!/usr/bin/env python

'''
Given a FASTA file with (multiple) sequences returns the transcript ids 
in which QRE pattern (1-2) is found.
1. Galarneau, A. & Richard, S. Target RNA motif and target mRNAs of the 
   Quaking STAR protein. Nat Struct Mol Biol 12, 691-698 (2005).
2. de Bruin, R. G. et al. Quaking promotes monocyte differentiation into 
   pro-atherogenic macrophages by controlling pre-mRNA splicing and gene 
   expression. Nat Commun 7, 10846 (2016).
@author: Enrico Gaffo
@email: enrico.gaffo@unipd.it
@license: any use of this script must be approved by Enrico Gaffo
'''

import re, argparse, gzip, sys
from Bio import SeqIO

def find_qre_pattern(sequence):
	'''
	Return the indexes and the stretch of the QRE pattern matches
	in a string readable format [start-end:stretch, start-end:stretch, ...]
	@param sequence: just a string where to search for the QRE pattern
	'''
	# store the QRE pattern for subsequent use
	qre_pattern = re.compile('[ACTG]ACTAA[CT][AGCT]{1,20}TAA[CT]')
	# use finditer in order to get also the positions of the match
	# with finditer() the string is scanned left-to-right, and 
	# matches are returned in the order found
	matches = qre_pattern.finditer(sequence)
	result = []
	# the pattern may be found multiple times within the same sequence
	for m in matches:
		# three element list
		matched_pos = [m.start(), m.end(), m.group(0)]
		result.append(matched_pos)
	# the result will be a list of three-element lists
	return result

if __name__ == '__main__':

	# just handle input parameters
	parser = argparse.ArgumentParser(description = \
		'''Given a FASTA file with (multiple) sequences returns the transcript ids '''
		'''in which QRE pattern (1) is found. '''
		'''[1]: de Bruin, R. G. et al. Quaking promotes monocyte differentiation into '''
		'''pro-atherogenic macrophages by controlling pre-mRNA splicing and gene '''
		'''expression. Nat Commun 7, 10846 (2016).''')
	parser.add_argument('infile', type = str, default = '-',
		 	    help = 'A gzipped FASTA file')
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
            # get sequence
            seq = str(trs.seq)
            geneID = trs.id
            # check if pattern is in current sequence
            matches = find_qre_pattern(seq)
            # if pattern matches one or more times in sequence
            if matches:
                matched_pos = []
        	for m in matches:
        	    qreStart = m[0]
        	    qreEnd   = m[1]
        	    qreString= m[2]
        	    matched_pos_str = str(qreStart)+ '-' +\
                                      str(qreEnd) + ':' +\
        	        	      str(qreString)
        	    matched_pos.append(matched_pos_str)
        
        	    # print results in TAB separated format
                print '\t'.join([geneID, str(len(matched_pos)), ';'.join(matched_pos)])
        
