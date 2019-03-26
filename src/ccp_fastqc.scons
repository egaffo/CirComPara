'''

This is a SConscript script that execute the FASTQC [1] program to analyze FASTA, 
FASTQ or BAM files.

Software dependencies:
 * FASTQC

When called from a SConscript it imports the following variables:
 * env
 * fastqc_readset
  
[1] http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
'''


import os, re

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_fastqc.Clone()
    #fastqc_readset
except NameError:
    vars = Variables('vars.py')
    vars.Add('READS', 'FASTQ read file to process, either in plain text or gzipped (.gz)', 
             'reads.fq')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)
    # These are the variables given from the command line when the SConscript is called
    # standalone
    fastqc_readset = File(env['READS'])


## fastqc_readset must be a Scons.File object
# FASTQC removes .gz, .fq, .fastq extensions (also recursively) from the input file name
# to name the output HTML and ZIP files
fastqc_readset = env['READSET']
readset_basename = re.sub('\.fq$|\.fastq$', '', re.sub('\.gz$', '', fastqc_readset.name))
fastqc_html_target = readset_basename + '_fastqc.html'
fastqc_zip_target = readset_basename + '_fastqc.zip'
## including the fastqc_data.txt file in the targets may result in 
## 'target doesn't exist' if empty fastq files are processed by FASTQC.
## We need to initialize the file.
fastqc_data_target = os.path.join(readset_basename + "_fastqc",
                                  "fastqc_data.txt")

fastqc_quality_cmd ='''echo "No reads in ${SOURCES[0]}" > ${TARGETS[0]} && '''\
                    '''echo "No reads in ${SOURCES[0]}" > ${TARGETS[4]} && '''\
                    '''fastqc $SOURCE -o ${TARGETS[0].dir} --extract > '''\
                    '''${TARGETS[2]} 2> ${TARGETS[3]} '''
fastqc_quality = env.Command([fastqc_html_target, 
                              fastqc_zip_target,
                              '${SOURCE.filebase}_fastqc.log', 
                              '${SOURCE.filebase}_fastqc.err', 
                              fastqc_data_target],
                             fastqc_readset, 
                             fastqc_quality_cmd)

## the following files comes from the unzipping of the zip'ed results
additional_files_to_clean = [os.path.join(readset_basename, 
                                          f) for f in ['fastqc.fo',
                                                       'fastqc_report.html',
                                                       'Icons', 
                                                       'Images', 
                                                       'summary.txt']]
Clean(fastqc_quality, additional_files_to_clean)

Return('fastqc_quality')
