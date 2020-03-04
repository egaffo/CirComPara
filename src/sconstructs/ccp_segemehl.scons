'''

This is a SConscript script that executes the tasks necessary to map
RNA-seq reads to a reference genome using 'segemehl' [1].


[1] Hoffmann, S. et al. 
    
    A multi-split mapping algorithm for circular RNA, splicing, 
    trans-splicing and fusion detection. 
    
    Genome Biology 15, R34 (2014).

Software dependencies:
 * Samtools-1.1

When called from a SConscript it imports the following variables:
 * env
 * segemehl_cpus
 * reads_to_map
 * sample_id
 * segemehl_genome_fasta
 * segemehl_index
 * segemehl_extra_params

'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_segemehl.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('SEGEMEHL_INDEX', 'The SEGEMEHL index', '')
    vars.Add('SEGEMEHL_PARAMS', 'SEGEMEHL extra parameters', '')
    vars.Add('GENOME_FASTA', '''The  path to genome. Point to folder with one '''
                             '''fasta file for each chromosome.''', '.')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    vars.Add('SAMPLE', 'Name of the sample', 'sample')
    
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

segemehl_mapping_dir = 'segemehl'

## MAP READS
#sam_sort_cmd = '''samtools view -u $(-@ $CPUS '''\
#               '''$) - | samtools sort $( -m $SAM_SORT_MM $) -@$( $CPUS '''\
#               '''$) -O 'sam' -T $( ${TARGETS[0].dir}_tmp_$SAMPLE '''\
#               '''$) | gzip > ${TARGETS[0].abspath} '''

## do not sort to prevent samtools failures to fail SAM output file
sam_sort_cmd = '''gzip > ${TARGETS[0].abspath} '''

segemehl_map_cmd = 'cd ${TARGETS[0].dir} && '
if len(env['READS']) > 1:
    segemehl_map_cmd = segemehl_map_cmd +\
                       '''segemehl.x $SEGEMEHL_PARAMS '''\
                       '''-i $SEGEMEHL_INDEX -d $GENOME_FASTA '''\
                       '''-q ${SOURCES[0].abspath} '''\
                       '''-p ${SOURCES[1].abspath} -S ''' +\
                       '''$(-t $CPUS $) '''\
                       '''-u ${TARGETS[1].filebase} | ''' + \
                       sam_sort_cmd
else:
    segemehl_map_cmd = segemehl_map_cmd +\
                       '''segemehl.x $SEGEMEHL_PARAMS '''\
                       '''-i $SEGEMEHL_INDEX -d $GENOME_FASTA '''\
                       '''-q ${SOURCES[0].abspath} -S '''\
                       '''$(-t $CPUS $) '''\
                       '''-u ${TARGETS[1].filebase} | ''' + \
                       sam_sort_cmd

segemehl_map_cmd = ' && '.join([segemehl_map_cmd, 
                                'cd ' + Dir('#').abspath, 
                                'gzip ${TARGETS[1].base}'])

segemehl_unmatched_target = "{}_unmatched.fastq.gz".format(env['SAMPLE'])

segemehl_map = env.Command([os.path.join(segemehl_mapping_dir, 
                                         "{}.unsorted.sam.gz".format(env['SAMPLE'])), 
                            os.path.join(segemehl_mapping_dir, 
                                         segemehl_unmatched_target),
                            os.path.join(segemehl_mapping_dir,
                                         '${SOURCES[0].filebase}.sngl.bed'),
                            os.path.join(segemehl_mapping_dir, 
                                         '${SOURCES[0].filebase}.mult.bed'),
                            os.path.join(segemehl_mapping_dir, 
                                         '${SOURCES[0].filebase}.trns.txt')],
                           env['READS'], 
                           segemehl_map_cmd)

# do not delete old alignment file untill the new one has been computed
env.Precious(segemehl_map[0]) 

## COUNT AND REPORT MAPPED READS
mapped_reads_target = os.path.join(segemehl_mapping_dir, 
                                   'segemehl_mapped_reads_count.txt')
mapped_reads_cmd    = '''zcat ${SOURCE} | samtools view -F 4 - '''\
                      '''| cut -f 1 | sort | uniq | wc -l > $TARGET'''
mapped_reads        = env.Command(mapped_reads_target, segemehl_map[0], mapped_reads_cmd) 

Clean('.', segemehl_mapping_dir)

results = {'ALIGNMENTS':      segemehl_map[0],
           'UNALIGNED_READS': segemehl_map[1],
           'SINGLE_SPLITS':   segemehl_map[2],
           'MULTI_SPLITS':    segemehl_map[3],
           'TRANS_SPLIT':     segemehl_map[4],
           'MAPPED_READS':    mapped_reads}

Return('results')
