'''

This is a SConscript script that executes the tasks necessary to map
RNA-seq reads to a reference genome using 'Hisat2' [1].

Kim D, Langmead B and Salzberg SL.
HISAT: a fast spliced aligner with low memory requirements.
Nature Methods 2015

Software dependencies:

When called from a SConscript it imports the following variables:
 * env
 * hisat2_cpus
 * reads_to_map
 * hisat2_index
 * hisat2_extra_params
 * sample_id

Returns:
    ['sample_hisat2.bam', 
     'sample_hisat2.log']
'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_hisat2.Clone()
    #CPUS = hisat2_cpus
    #READS = reads_to_map
    #HISAT2_INDEX = hisat2_index
    #HISAT2_PARAMS = hisat2_extra_params 
    #SAMPLE = sample_id

except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('HISAT2_INDEX', 'The Hisat2 index', '')
    vars.Add('HISAT2_PARAMS', 'Hisat2 extra parameters', '')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    vars.Add('SAMPLE', 'Name of the sample','sample')    

    cmdline_env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(cmdline_env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    # These are the variables given from the command line when the SConscript is called
    # standalone
    env = cmdline_env
    #CPUS = env['CPUS']
    #HISAT2_INDEX = env['HISAT2_INDEX'] # Assume you already have your Hisat2 genome index
    env['READS'] = env['READS'].split(',')
    #HISAT2_PARAMS = env['HISAT2_PARAMS']
    #SAMPLE = env['SAMPLE']
 
hisat2_out_dir = 'hisat2_out'
chdir_working_cmd = '''cd ''' + os.path.join(Dir('.').abspath, hisat2_out_dir)
hisat2_cmd = '''hisat2 -x $HISAT2_INDEX $HISAT2_PARAMS $(-p $CPUS $) '''

hisat2_targets_postfixes = [env['SAMPLE'] + '_hisat2.bam',
                            env['SAMPLE'] + '_hisat2.log']
if len(env['READS']) > 1:
    hisat2_cmd = hisat2_cmd + ''' -1 ${SOURCES[0].abspath} -2 ${SOURCES[1].abspath} '''
else:
    hisat2_cmd = hisat2_cmd + ''' -U $SOURCE.abspath --un-gz ${TARGETS[2].abspath} '''
    hisat2_targets_postfixes.append(env['SAMPLE'] + '_unmapped.fastq.gz')

chdir_parent_cmd = '''cd ''' + Dir('#').abspath

hisat2_targets = [os.path.join(Dir('.').abspath, 
                               hisat2_out_dir, f) for f in hisat2_targets_postfixes]

hisat2_cmd = hisat2_cmd + ''' 2> ${TARGETS[1].abspath} | samtools view -hu '''\
             '''-@ $( $CPUS $) - | samtools sort $( -m $SAM_SORT_MM $) -O 'bam' -@ $( $CPUS $) '''\
             '''-T $(hisat2_''' + env['SAMPLE'] + '''$) > ${TARGETS[0].abspath}'''

hisat2 = env.Command(hisat2_targets, 
                     env['READS'], 
                     ' && '.join([chdir_working_cmd, 
                                  hisat2_cmd,
                                  chdir_parent_cmd]))

env.Precious(hisat2[0]) # do not delete old alignment file until the new one has been computed

## make BAM index (useful for downstream analysis tools)
bam_index_target = '''${SOURCE.abspath}.bai'''
bam_index_cmd = '''samtools index $SOURCE ${SOURCE.path}.bai'''
bam_index = env.Command(bam_index_target, hisat2[0], bam_index_cmd)

Clean('.', hisat2_out_dir)

results = {'BAM' : hisat2[0],
           'LOG' : hisat2[1],
           'BAM_INDEX' : bam_index,
           'UNMAPPED_READS': None}

if len(env['READS']) < 2:
    results['UNMAPPED_READS'] = hisat2[2]


Return('results')
