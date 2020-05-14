'''

This is a SConscript script that executes the tasks necessary to find 
circular RNA junctions from RNA-seq data, according to the procedure 
described in:

    Gao, Y., Wang, J. & Zhao, F. 
    
    CIRI: an efficient and unbiased algorithm for de novo circular RNA
    identification. 
    
    Genome Biology 16, 4 (2015).

Download the scripts from http://sourceforge.net/projects/ciri

Returns: 
 [sampleid_ciri.out]

'''

import os


Import('*')

try:
    # these variables can be passed with 'exports' when calling this SConscript
    # from another SConscript
    env = env_ciri.Clone()

except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('GENOME_FASTA', '''The  path to genome. Point to folder with one'''
                             ''' fasta file for each chromosome.''', '.')
    vars.Add('SAMPLE', 'Name of the sample', '')
    vars.Add('CIRI', 'The full path to the CIRI_vx.x.pl perl script', '')
    vars.Add('ANNOTATION', 'The full path to the GTF/GFF annotation file', '')
    vars.Add('CIRI_EXTRA_PARAMS', 'CIRI additional parameters', '')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

out_dir = 'ciri_out'
chdir_working_cmd  = 'cd ' + os.path.join(Dir('.').abspath, out_dir)
chdir_previous_cmd = 'cd ' + Dir('#').abspath

#ciri_parameters = env['CIRI_EXTRA_PARAMS']
if isinstance(env['CIRI_EXTRA_PARAMS'], basestring):
    env['CIRI_EXTRA_PARAMS'] = env['CIRI_EXTRA_PARAMS'].split()

if env['ANNOTATION']:
    #ciri_parameters = ciri_parameters + ' -A ' + env['ANNOTATION']
    env['CIRI_EXTRA_PARAMS'].extend(['-A', env['ANNOTATION']])

ciri_target = os.path.join(out_dir, env['SAMPLE'] + '_ciri.out')

## Run CIRI on mappings
## multithread -T option works from v2 of CIRI
ciri_cmd = ' && '.join(['zcat ${SOURCES[0].abspath} > ${SOURCES[0].filebase}.temp',
                        'perl $CIRI -T $( $CPUS $) -I ${SOURCES[0].filebase}.temp '\
                        '-O ${TARGETS[0].file} -F $GENOME_FASTA $CIRI_EXTRA_PARAMS',
                        'rm ${SOURCES[0].filebase}.temp'])

## Run the commands
ciri = env.Command([ciri_target],
                   env['BWA_ALIGN'], 
                   ' && '.join([chdir_working_cmd, ciri_cmd, chdir_previous_cmd]))

## Collect backsplice read names
#bks_reads_cmd = '''cut -f 12 ${SOURCES[0]} | '''\
#                '''grep -v "junction_reads_ID" | '''\
#                '''sed -E "s/,$//" | sed -E "s/,/\\n/g" | '''\
#                '''sort | uniq -c | '''\
#                '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
#                '''sort -k1,1nr > ${TARGETS[0]} '''
#
#bks_reads = env.Command([os.path.join(out_dir, 
#                                      env['SAMPLE'] + '.ciri.bks.reads')], 
#                        [ciri], 
#                        bks_reads_cmd)

bks_reads_tagets = [os.path.join(out_dir, f) for f in 
                                ['${SAMPLE}.circular.reads.bed.gz', 
                                 '${SAMPLE}.ciri.bks.reads']]
bks_reads_cmd = 'get_ciri_bks_reads.R -i ${SOURCES[0]} -b ${TARGETS[0]} -l ${TARGETS[1]}'
bks_reads = env.Command(bks_reads_tagets, ciri, bks_reads_cmd)

## Clean log files, etc.
Clean(ciri, [ciri_target + '.log', os.path.join(Dir('.').abspath, 'CIRIerror.log')])

results = {'CIRCRNAS':    ciri[0],
           #'CIRC_SN_BED': bed,
           'BKS_READS':   bks_reads[1],
           'CIRC_READS':  bks_reads[0]}

Return('results')
