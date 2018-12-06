#!/usr/bin/env python

import csv, argparse, sys

def format_gtf_line(chrom, source, feature, start, end, score, strand, frame, 
                    gene_id, transcript_id, extra_field = ''):
    group = 'gene_id "' + gene_id + '"; transcript_id "' +\
            transcript_id + '"; ' + extra_field
    line = '\t'.join([chrom, source, feature, start, end, score, strand, frame, group])
    return line

def format_circexplorer(line, outformat, ce_version):
    fields = line 
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        start   = str(int(float(fields[2])) + 1) # CIRCexplorer gives BED12 coordinates
        end     = fields[3]
        strand  = fields[6]
        score   = fields[13]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, ce_version, 'backsplice', 
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
        start   = str(int(float(fields[2])) + 1) # findcirc gives BED coordinates
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

def format_testrealign(line, outformat):
    fields = line
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        start   = str(int(float(fields[2])) + 1) # testrealign gives BED coordinates
        end     = fields[3]
        strand  = fields[5]
        score   = fields[6]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'testrealign', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')

    elif outformat == '':
        outline = 'BED6 formatting not yet immplemented'
    return outline

def format_dcc(line, outformat):
    fields = line
    if outformat == 'gtf':
        sample = fields[0]
        chrom   = fields[1]
        #start   = str(int(float(fields[2])) + 1) # DCC gives BED coordinates
        # DCC gives BED coordinates, but seems it does not respect the 0 based position counting
        start   = str(int(float(fields[2]))) 
        end     = fields[3]
        strand  = fields[4]
        score   = fields[5]
        gene_id = chrom + ':' + start + '-' + end + ':' + strand
        transcript_id = gene_id + '.' + sample
        
        outline = format_gtf_line(chrom, 'dcc', 'backsplice', 
                                  start, end, score, strand, '.', 
                                  gene_id, transcript_id, 'sample_id "' + sample + '";')

    elif outformat == 'bed6':
        outline = 'BED6 formatting not yet immplemented'
    return outline

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = '''Convert Junk2-circpipe circRNA '''\
                                     '''result collection table into an annotation file''', 
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', default = '-', help = 'A csv file or piped stream (default -)')
    parser.add_argument('-p', '--program', type = str, 
                        choices = ['ciri', 'circexplorer', 'findcirc', 'testrealign', 
                                   'circexplorer2_star', 
                                   'circexplorer2_bwa', 
                                   'circexplorer2_segemehl',
                                   'circexplorer2_tophat_pe',
                                   'circexplorer2_tophat',
                                   'circexplorer2_mapsplice',
                                   'dcc'], 
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
        
        if args.program in ['circexplorer', 'circexplorer2_star',
                            'circexplorer2_bwa', 'circexplorer2_segemehl',
                            'circexplorer2_tophat_pe', 'circexplorer2_tophat', 
                            'circexplorer2_mapsplice']:
            outline = format_circexplorer(line, args.outformat, args.program)
        elif args.program == 'ciri':
            outline = format_ciri(line, args.outformat)
        elif args.program == 'findcirc':
            outline = format_findcirc(line, args.outformat)
        elif args.program == 'testrealign':
            outline = format_testrealign(line, args.outformat)
        elif args.program == 'dcc':
            outline = format_dcc(line, args.outformat)

        print outline

