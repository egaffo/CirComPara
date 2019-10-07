#!/usr/bin/env python

import argparse, sys


if __name__ == '__main__':
    
    desc = ''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('infile', 
                        default = '-', 
                        help = 'Input GTF file. Set - for stdin stream.')
    parser.add_argument('-o', '--outfile', 
                        type = str, 
                        default = '-',
                        help = 'Output filename. Default to stdout.')
    parser.add_argument('-t', '--trim', 
                        action = 'store_true',
                        help = 'Trim the GTF attribute field and keep only '\
                               'gene_id.')

    args = parser.parse_args()

    if args.infile == '-':
        infile = sys.stdin
    else:
        infile = open(args.infile, 'r')

    if args.outfile == '-':
        outfile = sys.stdout
    else:
        outfile = open(args.outfile, 'w')

    for inline in infile:
        outline = inline.split('\t')
        attribute = outline[8]

        if args.trim:
            for field in outline[8].split(';'):
                if 'gene_id' in field.split():
                    attribute = field + ';\n'

        start_out = '\t'.join(outline[0:2] +
                             ['start'] +
                             [outline[3]] +
                             [outline[3]] +
                             outline[5:8] +
                             [attribute])
        stop_out  = '\t'.join(outline[0:2] +
                              ['stop'] +
                              [outline[4]] +
                              [outline[4]] +
                              outline[5:8] +
                              [attribute])
        outfile.write(start_out)
        outfile.write(stop_out)

    infile.close()
    outfile.close()
        
        
