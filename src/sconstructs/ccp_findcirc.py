'''

This is a SConscript script that executes the tasks necessary to find 
circular RNA junctions from RNA-seq data, according to the procedure 
described in:

  Nature. 2013 Mar 21;495(7441):333-8. 
  doi: 10.1038/nature11928. Epub 2013 Feb 27.
      
  Circular RNAs are a large class of animal RNAs with regulatory potency.
          
  Memczak S, Jens M, Elefsinioti A, Torti F, Krueger J, Rybak A, Maier L,
  Mackowiak SD, Gregersen LH, Munschauer M, Loewer A, Ziebold U, 
  Landthaler M, Kocks C, le Noble F, Rajewsky N.

Download the scripts from http://www.circbase.org/download/find_circ.tar.gz

Software dependencies: 
 * Bowtie2-2.2.4
 * Pysam
 * find_circ.py v1.2

When called from a SConscript it imports the following variables:
 * CPUS
 * BOWTIE2_INDEX
 * GENOME_FASTA
 * READS
 * SAMPLE
 * FINDCIRC_EXTRA_PARAMS

'''
import os

Import('*')

try:
    # these variables can be passed with 'exports' when calling this SConscript
    # from another SConscript
    env = env.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('BOWTIE2_INDEX', 'The Bowtie2 index', 'bt2_cdr1as_locus')
    vars.Add('GENOME_FASTA', 'The  path to genome. Point to folder with one fasta file for each chromosome.', '.')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    vars.Add('SAMPLE', 'Name of the sample', 'cdr1as_test_')
    vars.Add('FINDCIRC_EXTRA_PARAMS', 'Extra parameters to be passed to find_circ.py script.'\
             'F.i. --stranded --strandpref --halfunique --noncanonical', 
             '')
    vars.Add('BOWTIE2_PARAMS', 'Extra parameters to pass to Bowtie2 in addition to'\
            '-p $CPUS --reorder --score-min=C,-15,0 -q', '')
    env = Environment(ENV=os.environ,
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

out_dir = 'find_circ_out'

env.SetDefault(FINDCIRC_EXTRA_PARAMS = '')

best_qual_val = 40
if '--best-qual' in env['FINDCIRC_EXTRA_PARAMS']:
    best_qual_idx = env['FINDCIRC_EXTRA_PARAMS'].index('--best-qual') + 1
    best_qual_val = env['FINDCIRC_EXTRA_PARAMS'].pop(best_qual_idx)
    env['FINDCIRC_EXTRA_PARAMS'].remove('--best-qual')

## Findcirc considers only single-end reads. If CirComPara paired-end
## mode is enabled, we need to make consistent the direction of the 
## read mates before concatenating them. Moreover, we need to make 
## consistent also the read names for post-processing method comparison
## by fixing the modified read names
if env['CIRC_PE_MAPPING']:
    revcomp_cmd = '''fastq_rev_comp.py -f $SOURCE | '''\
                  '''sed 's/revcomp_of_//' | '''\
                  '''gzip -c > $TARGET'''
    revcomp = env.Command(os.path.join(out_dir, 
                                       env['SAMPLE'] + 
                                       '_1.revcomp.fq.gz'), 
                          env['READS'][0], 
                          revcomp_cmd) 
    env['READS'][0] = revcomp[0]
    
cat_cmd = 'cat'
if File(env['READS'][0]).path.endswith('.gz'):
    cat_cmd = 'zcat'

find_circ_cmd = 'unmapped2anchors.py -Q <( ' + cat_cmd + ' ${SOURCES} ) |' + \
                ' bowtie2 $BOWTIE2_PARAMS $( -p $CPUS $) --reorder --score-min=C,-15,0 -q -x ' + \
                ' $BOWTIE2_INDEX -U - 2> ${TARGETS[2]} | find_circ.py ' + \
                ' -G $GENOME_FASTA -p ${SAMPLE}_ -s ${TARGETS[3]}' + \
                ' -R ${TARGETS[1]} $FINDCIRC_EXTRA_PARAMS > ${TARGETS[0]}'


find_circ = env.Command([os.path.join(out_dir, 'sites.bed'), 
                         os.path.join(out_dir, 'sites.reads'),
                         os.path.join(out_dir, 'bt2_secondpass.log'),
                         os.path.join(out_dir, 'find_circ.log')], 
                        env['READS'], 
                        find_circ_cmd)

# Select backsplices just by their tag. Additional filters, like minimum read number
# and length are not performed
# NB: the commands within curly braces are to prevent grep to return error status
# when no match and/or no input
#filter_circ_cmd = "grep CIRCULAR ${SOURCES[0]} | " +\
#                  " { grep UNAMBIGUOUS_BP || true; } | "+\
#                  " { grep ANCHOR_UNIQUE || true; } > ${TARGET}"
filter_circ_cmd = 'filter_findcirc_res.R -i $SOURCE -o $TARGET -q ' +\
                  str(best_qual_val)
filter_circ = env.Command([os.path.join(out_dir, 'circ_candidates.bed')], 
                          find_circ[0], 
                          filter_circ_cmd)

Clean('.', out_dir)

## Collect backsplice read names
## use the sites.reads file: consider only the FASTA headers
## in which it is reported the read name as a second field
## with attached an underscore and a number
## N.B.: count only circRNAs in circ_candidates.bed
#bks_reads_cmd = '''grep ">" ${SOURCES[0]} | '''\
#                '''grep -f <(cut -f 4 ${SOURCES[1]}) | '''\
#                '''sed -E "s/>[^ ]+ ([^_]+)_[0-9]+/\\1/g" | '''\
#                '''sort | uniq -c | '''\
#                '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
#                '''sort -k1,1nr > ${TARGETS[0]} '''
#
#bks_reads = env.Command([os.path.join(out_dir, env['SAMPLE'] + '.fc.bks.reads')], 
#                        [find_circ[1], filter_circ], 
#                        bks_reads_cmd)

bks_reads_tagets = [os.path.join(out_dir, f) for f in 
                                ['${SAMPLE}.circular.reads.bed.gz', 
                                 '${SAMPLE}.fc.bks.reads']]
bks_reads_cmd = 'get_findcirc_bks_reads.R -s ${SOURCES[0]} '\
                '-c ${SOURCES[1]} -r ${TARGETS[0]} -l ${TARGETS[1]}'
bks_reads = env.Command(bks_reads_tagets, 
                        [find_circ[1], filter_circ[0]], 
                        bks_reads_cmd)

results = {'CIRCRNAS':    filter_circ[0],
           #'CIRC_SN_BED': bed,
           'BKS_READS':   bks_reads[1],
           'CIRC_READS':  bks_reads[0]}

#Return('find_circ filter_circ mapped_reads bks_reads')
Return('results')
