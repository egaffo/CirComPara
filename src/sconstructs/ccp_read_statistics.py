'''

This is a SConscript script that performs various scripts to analyse FASTQ files.
It relies on FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) and custom
R scripts.
Moreover, it executes quality cehck with FASTQC

Software dependencies:
 * FASTX-Toolkit v0.0.14
 * FASTQC
 * CIRCOMPARA custom R scripts

When called from a SConscript it imports the following variables:
 * env
 * read_statistics_readset
  
'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_read_statistics.Clone()
    readset = env['READSET']#read_statistics_readset
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('READS', 'FASTQ read file to process, either in plain text or gzipped (.gz)', 
             'reads.fq')
    vars.Add('READSTAT_METHODS', 'Comma separated list of methods to use for read statistics. '\
             'Currently supported: fastqc,fastx', 'fastqc')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)
    # These are the variables given from the command line when the SConscript is called
    # standalone
    readset = env['READS']

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src')

env.SetDefault(READSTAT_METHODS = 'fastqc')

statistics_dir = 'fastx_stats'

readset_basename = os.path.splitext(os.path.basename(readset))[0]
readset_ext = os.path.splitext(os.path.basename(readset))[1]

if readset_ext == '.gz':
    pre_cmd = 'zcat $SOURCE | '
else:
    pre_cmd = 'cat $SOURCE | '

stats = []
## COMPUTE STATISTICS
if('fastx' in env['READSTAT_METHODS'].split(',')):

    fastx_quality_stats_cmd = 'fastx_quality_stats -N > $TARGET'
    fastx_quality_stats = env.Command(os.path.join(statistics_dir, 
                                                  readset_basename + '.fastx_quality_stats.tsv'), 
                                     readset, pre_cmd + fastx_quality_stats_cmd)
    stats.append(fastx_quality_stats)
    
    ## DRAW PLOTS:
    fastx_quality_plot_cmd = 'fastx_stats_make_boxplot.R -i $SOURCE -o $TARGET'
    fastx_quality_plot= env.Command(os.path.join(statistics_dir,
                                                 readset_basename + '.fastx_quality_boxplot.svg'),
                                     fastx_quality_stats, fastx_quality_plot_cmd)
    stats.append(fastx_quality_plot)
    
    length_plot_cmd = 'fastx_stats_make_length_plot.R -i $SOURCE -o $TARGET'
    length_plot = env.Command(os.path.join(statistics_dir, readset_basename + '.length_plot.svg'),
                              fastx_quality_stats, length_plot_cmd)
    stats.append(length_plot)
    
    nucleotide_plot_cmd = 'fastx_stats_make_nucleotide_plot.R -i $SOURCE -o $TARGET'
    nucleotide_plot = env.Command(os.path.join(statistics_dir, 
                                               readset_basename + '.nucleotide_plot.svg'),
                                  fastx_quality_stats, nucleotide_plot_cmd)
    stats.append(nucleotide_plot)
    
    Clean('.', statistics_dir)

ccp_fastqc = []
if('fastqc' in env['READSTAT_METHODS'].split(',')):
    ## perform FASTQC analyses (default)
    fastqc_dir   = 'fastqc_stats'
    env_fastqc = env.Clone()
    #fastqc_readset = File(readset)
    env_fastqc['READSET'] = File(readset)
    ccp_fastqc = env.SConscript(os.path.join(fastqc_dir, 'ccp_fastqc.py'), 
                                  src_dir = SRC_DIR,
                                  variant_dir = fastqc_dir, duplicate = 0,
                                  exports = 'env_fastqc')
    Clean('.', fastqc_dir)

Return('stats ccp_fastqc')
