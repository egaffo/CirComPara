#!/usr/bin/env python

import csv, argparse, sys, re

def fix_coords_and_name(fields, bed, name, source_name = 'junk2'):
    ## convert BED coordinates (0-based and end position not included) 
    ## into GTF (1-based and end position included). NB: 'end' values remain the same!
    bed2gtf = [[bed[0][0], int(bed[0][1]) +1, int(bed[0][2])], 
               [bed[1][0], int(bed[1][1]) +1, int(bed[1][2])]]
    if not name:
        name = 'translocated_' + '_'.join([str(p) for p in bed2gtf[0][0:2]]) + '_' +\
               '_'.join([str(p) for p in bed2gtf[1][0:2]])
    tr_offset = bed2gtf[0][2] - bed2gtf[0][1] + 1 # the lenght of left translocated chromosome
    chrom   = fields[0]
    source  = fields[1]
    feature = fields[2]
    start   = fields[3]
    end     = fields[4]
    strand  = fields[6]
    group   = fields[8]
    if chrom == bed2gtf[0][0] or chrom == bed2gtf[1][0]:
        group   = re.sub(r'transcript_id "', 'transcript_id "'+ name + '_', group)
        group   = re.sub(r'gene_id "', 'gene_id "'+ name + '_', group)
    if chrom == bed2gtf[0][0]:
        # remove translocation offset from coordinates
        start   = str(max(1, (int(fields[3]) - bed2gtf[0][1]) + 1))
        end     = str(min((int(fields[4]) - bed2gtf[0][1]) + 1, bed2gtf[0][1]))
        #strand  = '+' # is the strand preserved after translocation? Leave commented if 'yes'.
        chrom   = name
        source  = source_name
    elif chrom == bed2gtf[1][0]:
        # remove translocation offset from coordinates and add offset of the prepended chromosome
        start   = str(max(1, (int(fields[3]) - bed2gtf[1][1]) + 1) + tr_offset)
        end     = str(min((int(fields[4]) - bed2gtf[1][1]) + 1, bed2gtf[1][1]) + tr_offset)
        #strand  = '+' 
        chrom   = name
        source  = source_name
    score   = fields[5]
    frame   = fields[7]
    outline = '\t'.join([chrom, source, feature, start, end, score, strand, frame, group])
    return outline


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = '''Adjust annotation coordinates to be '''\
                                     '''relative to the translocated chromosome.''', 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', default = '-', help = 'A GTF file or piped stream (default -)')
    parser.add_argument('-b', '--bed', required = True, dest = 'bed', 
                        help = 'The BED file with translocation coordinates')
    parser.add_argument('-n', '--name', type = str, dest = 'chrname', default = None,
                        help = 'The name of the translocated chromosome')
    parser.add_argument('-s', '--source_name', type = str, dest = 'srcname', default = 'junk2',
                        help = 'The name of the program that is modifying the annotation')
    args = parser.parse_args()

    # will result in a list of two 3-elements-lists representing bed coordinates
    bed = [l for l in csv.reader(open(args.bed, 'r'), delimiter='\t')]

    if args.input == '-':
        f = sys.stdin
    else:
        f = open(args.input, 'r')

    linefields = csv.reader(f, delimiter='\t')
    for line in linefields:
        if not line[0].startswith('#'):
            outline = fix_coords_and_name(line, bed, args.chrname, args.srcname)
            print outline
        else:
            print '\t'.join(line)


