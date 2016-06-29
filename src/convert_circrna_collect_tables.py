#!/usr/bin/env python

import csv, argparse, sys

def format_gtf_line(chrom, source, feature, start, end, score, strand, frame, 
                    gene_id, transcript_id, extra_field = ''):
    group = 'gene_id "' + gene_id + '"; transcript_id "' +\
            transcript_id + '"; ' + extra_field
    line = '\t'.join([chrom, source, feature, start, end, score, strand, frame, group])
    return line

def format_circexplorer(line, outformat):
    fields = line 
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        start   = str(int(fields[2]) + 1) # CIRCexplorer gives BED12 coordinates
        end     = fields[3]
        strand  = fields[6]
        score   = fields[13]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'circexplorer', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')
    elif outformat == 'bed6':
        outline = 'BED6 formatting not yet immplemented'
    return outline

def format_ciri(line, outformat):
    fields = line
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[2]
        start   = fields[3] # CIRI gives GFF coordinates
        end     = fields[4]
        strand  = fields[11]
        score   = fields[5]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'ciri', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')
    elif outformat == 'bed6':
        outline = 'BED6 formatting not yet immplemented'
    return outline

def format_findcirc(line, outformat):
    fields = line
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        start   = str(int(fields[2]) + 1) # findcirc gives BED coordinates
        end     = fields[3]
        strand  = fields[6]
        score   = fields[5]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'findcirc', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')

    elif outformat == 'bed6':
        outline = 'BED6 formatting not yet immplemented'
    return outline

def format_segecirc(line, outformat):
    fields = line
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        start   = str(int(fields[2]) + 1) # segecirc gives BED coordinates
        end     = fields[3]
        strand  = fields[5]
        score   = fields[6]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'segecirc', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')

    elif outformat == '':
        outline = 'BED6 formatting not yet immplemented'
    return outline


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = '''Convert Junk2-circpipe circRNA '''\
                                     '''result collection table into an annotation file''', 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', default = '-', help = 'A csv file or piped stream (default -)')
    parser.add_argument('-p', '--program', type = str, 
                        choices = ['ciri', 'circexplorer', 'findcirc', 'segecirc'], 
                        required = True, dest = 'program', 
                        help = 'The program that generated the input file')
    parser.add_argument('-f', '--format', type = str, 
                        choices = ['gtf', 'bed6'], default = 'gtf', dest = 'outformat', 
                        help = 'The format to convert into')
    args = parser.parse_args()

    if args.input == '-':
        f = sys.stdin
    else:
        f = open(args.input, 'r')
    
    next(f) # skip header line

    linefields = csv.reader(f, delimiter=',', quotechar='"')
    for line in linefields:
  
        outline = ''
        
        if args.program == 'circexplorer':
            outline = format_circexplorer(line, args.outformat)
        elif args.program == 'ciri':
            outline = format_ciri(line, args.outformat)
        elif args.program == 'findcirc':
            outline = format_findcirc(line, args.outformat)
        elif args.program == 'segecirc':
            outline = format_segecirc(line, args.outformat)

        print outline

