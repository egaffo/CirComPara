'''
This script processes RNA-seq reads by means of Trimmomatic:

    Bolger, A. M., Lohse, M. & Usadel, B. 
    Trimmomatic: a flexible trimmer for Illumina sequence data. 
    Bioinformatics 30, 2114-2120 (2014).

    http://www.usadellab.org/cms/?page=trimmomatic

Required variables to export when calling this script from a SConscript:
 * trimmomatic_cpus
 * trimmomatic_adapter_file
 * trimmomatic_extra_params
 * trimmomatic_reads

Note that you must set in your environment the TRIMMOMATIC variable that specifies
the full path to the Trimmomatic JAR file.
'''

import os, itertools

Import('*')

try:
    env                 = env
    CPUS                = trimmomatic_cpus
    ADAPTER_FILE        = trimmomatic_adapter_file
    TRIMMOMATIC_PARAMS  = trimmomatic_extra_params
    READS               = trimmomatic_reads
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', '''Number of parallel jobs to execute''', '4')
    vars.Add('ADAPTER_FILE', 
             '''The file with 3' adapter. Leave empty if no adapter trimming needed''', '')
    vars.Add('TRIMMOMATIC_PARAMS', 
             '''Trimmomatic parameters. Note that you must declare the steps to perform '''
             '''(e.g. in default value)''', 
             '''MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30''')
    vars.Add('READS', '''Reads to process. Gzipped (.gz) files accepted.'''
                      ''' Comma separated list for paired-end''', 'reads.fa')

    env = Environment(ENV=os.environ, variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    # These are the variables given from the command line when the SConscript is called
    # standalone
    CPUS                = env['CPUS']
    ADAPTER_FILE        = env['ADAPTER_FILE']
    TRIMMOMATIC_PARAMS  = env['TRIMMOMATIC_PARAMS']
    READS               = env['READS'].split(',')

trimmomatic_dir = 'trimmomatic'

if not ADAPTER_FILE.strip() == '':
    adapter_removal = ' ILLUMINACLIP:' + ADAPTER_FILE + ':2:30:10 '
else:
    adapter_removal = ''
if not TRIMMOMATIC_PARAMS=='':
    trimmomatic_steps = TRIMMOMATIC_PARAMS
else:
    trimmomatic_steps = 'MAXINFO:40:0.5 LEADING:20 TRAILING:20 '\
                        'SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30'

TRIMMOMATIC_JAR = env['ENV']['TRIMMOMATIC']

reads = READS
if len(reads)>1:
    ## MANAGE PAIRED END READS
    read1_filebase = os.path.splitext(os.path.basename(reads[0].replace('.gz', '')))[0]
    read2_filebase = os.path.splitext(os.path.basename(reads[1].replace('.gz', '')))[0]
    targets = [os.path.join(trimmomatic_dir, '{}.fq.P.qtrim.gz'.format(read1_filebase)),
               os.path.join(trimmomatic_dir, '{}.fq.U.qtrim.gz'.format(read1_filebase)),
               os.path.join(trimmomatic_dir, '{}.fq.P.qtrim.gz'.format(read2_filebase)),
               os.path.join(trimmomatic_dir, '{}.fq.U.qtrim.gz'.format(read2_filebase)),
               os.path.join(trimmomatic_dir, 'trimmomatic.log')]
    trimmomatic_cmd = '$(java -jar ' + TRIMMOMATIC_JAR + '$) PE $(-threads ' +\
                      CPUS + '$) ${SOURCES[0]} ${SOURCES[1]} '\
                      '${TARGETS[0]} ${TARGETS[1]} ${TARGETS[2]} ${TARGETS[3]} ' +\
                      adapter_removal + trimmomatic_steps + ' 2> ${TARGETS[4]}'
else:
    ## MANAGE SINGLE END READS
    read1_filebase = os.path.splitext(os.path.basename(reads[0].replace('.gz', '')))[0]
    targets = [os.path.join(trimmomatic_dir, '{}.fq.SE.fq.gz'.format(read1_filebase)),
               os.path.join(trimmomatic_dir, 'trimmomatic.log')]
    trimmomatic_cmd = '$(java -jar ' + TRIMMOMATIC_JAR + '$) SE $(-threads ' +\
                      CPUS + '$) ${SOURCES[0]} ${TARGETS[0]} ' + \
                      adapter_removal + trimmomatic_steps + ' 2> ${TARGETS[1]}'

clean_reads = env.Command(targets, reads, trimmomatic_cmd)

Clean('.', trimmomatic_dir)
Return('clean_reads')
