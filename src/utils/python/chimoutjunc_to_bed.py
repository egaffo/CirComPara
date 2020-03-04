#!/usr/bin/env python

import os, sys, argparse, gzip

if __name__ == '__main__':

    usage = """
           Convert STAR's Chimeric.out.junction into BED coordinates. Output in
           stdout.
            """

    parser = argparse.ArgumentParser(description = usage, 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infilename', default = '-', required = True,
                        help = 'Chimeric.out.junction file path')
    args = parser.parse_args()


    if args.infilename == '-':
        junc_f = sys.stdin
    else:
        junc_f = open(args.infilename, 'r')
    
    for line in junc_f:
        flag = int(line.split()[6])
        if flag < 0:
            continue
        chr1, site1, strand1, chr2, site2, strand2, score1, score2, score3, read_id = line.split()[:10]
        if chr1 != chr2 or strand1 != strand2:
            continue
        if strand1 == '+':
            start = int(site2)
            end   = int(site1) - 1
        else:
            start = int(site1)
            end   = int(site2) - 1
        if start > end:
            continue
        else:
            bed_line = '\t'.join([chr1, str(start), str(end), read_id, score1, strand1])
            print bed_line

