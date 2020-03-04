'''

This is a SConscript script that executes the tasks necessary to map
RNA-seq reads to a reference genome using 'TopHat2' [1].

[1] 1.Kim, D. et al. 
    TopHat2: accurate alignment of transcriptomes in the presence 
    of insertions, deletions and gene fusions. 
    Genome Biology 14, R36 (2013).

Returns:
{ 'ALIGNMENTS' : 'accepted_hits.bam',
  'JUNCTIONS' : 'junctions.bed',
  'INSERTIONS': 'insertions.bed',
  'DELETIONS' : 'deletions.bed',,
  'UNMAPPED'  : 'unmapped.bam',
  'MAPPED_READS': 'tophat_mapped_reads_count.txt'
}

'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_tophat.Clone()

except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('TOPHAT_INDEX', 'The TopHat index', '')
    vars.Add('TOPHAT_PARAMS', 'TopHat extra parameters', '')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    
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
    env.Replace(READS = env['READS'].split(','))
 
tophat_out_dir = 'tophat_out'
tophat_cmd = '''tophat2 -o ${TARGETS[0].dir}''' \
             ''' $TOPHAT_PARAMS $(-p $CPUS $) $TOPHAT_INDEX '''

READS = [File(f) for f in env['READS']]

if len(READS) > 1:
    tophat_cmd = tophat_cmd + ''' ${SOURCES[0].abspath} ${SOURCES[1].abspath} '''
else:
    tophat_cmd = tophat_cmd + ''' $SOURCE.abspath '''

tophat_targets = [os.path.join(tophat_out_dir, f) for f in ['accepted_hits.bam', 
                                                            'junctions.bed',
                                                            'insertions.bed',
                                                            'deletions.bed',
                                                            'unmapped.bam']]

tophat_src = READS
## update source file dependencies for TopHat:
## annotation (file content) may change
if '--GTF' in env['TOPHAT_PARAMS']:
    tophat_src.append(File(env['TOPHAT_PARAMS'][env['TOPHAT_PARAMS'].index('--GTF')+1]))

tophat = env.Command(tophat_targets, 
                     tophat_src, 
                     tophat_cmd)


# do not delete old alignment file untill the new one has been computed
env.Precious(tophat[0])

results = { 'ALIGNMENTS' : tophat[0],
            'JUNCTIONS' : tophat[1],
            'INSERTIONS': tophat[2],
            'DELETIONS' : tophat[3],
            'UNMAPPED'  : tophat[4]
            }

## COUNT AND REPORT MAPPED READS
mapped_reads_target = os.path.join(tophat_out_dir, 'tophat_mapped_reads_count.txt')
mapped_reads_cmd    = '''samtools view -F 4 $SOURCE '''\
                      '''| cut -f 1 | sort | uniq | wc -l > $TARGET'''
mapped_reads        = env.Command(mapped_reads_target, results['ALIGNMENTS'], mapped_reads_cmd)

results['MAPPED_READS'] = mapped_reads

Clean('.', tophat_out_dir)

Return('results')
