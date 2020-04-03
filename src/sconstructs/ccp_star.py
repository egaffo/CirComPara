'''

This is a SConscript script that executes the tasks necessary to map
RNA-seq reads to a reference genome using 'Star' [1].

[1] Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., and Gingeras, T.R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21.

Software dependencies:
 * 

When called from a SConscript it imports the following variables:
 * env_star

Returns:
 [['Aligned.sortedByCoord.out.bam', 
   'Chimeric.out.junction',
   'Log.out', 
   'Log.progress.out',
   'Log.final.out', 
   'SJ.out.tab'],
  'STAR_mapped_reads_count.txt']

'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = star_env.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('STAR_INDEX', 'The Star index directory', '/path/to/index/dir')
    vars.Add('STAR_PARAMS', 'TopHat extra parameters', '')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    vars.Add('ANNOTATION', 'Genome annotation GTF file', 'genes.gtf')
    
    env = Environment(ENV = os.environ,
                      variables = vars)
    
    Help(vars.GenerateHelpText(cmdline_env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

star_out_dir = 'star_out'
chdir_working_cmd = '''cd ''' + os.path.join(Dir('.').abspath, star_out_dir) 
star_cmd = '''STAR --outSAMstrandField intronMotif '''\
            '''--outSAMtype BAM SortedByCoordinate '''\
            '''$STAR_PARAMS $(--runThreadN $CPUS $) '''\
            '''--genomeDir $STAR_INDEX'''

READS = [File(f) for f in env['READS']]

star_src = READS

if len(READS)>1:
    star_cmd = star_cmd + ''' --readFilesIn ${SOURCES[0].abspath} ${SOURCES[1].abspath} '''
else:
    star_cmd = star_cmd + ''' --readFilesIn $SOURCE.abspath '''

if READS[0].abspath.endswith('.gz'):
    star_cmd = star_cmd + ''' --readFilesCommand zcat '''

if not env['ANNOTATION'] == '':
    star_cmd = star_cmd + ''' --sjdbGTFfile $ANNOTATION'''
    star_src.append(env['ANNOTATION'])

chdir_parent_cmd = '''cd ''' + Dir('#').abspath

star_out_files = ['Aligned.sortedByCoord.out.bam', 
                  'Chimeric.out.junction',
                  'Log.out', 'Log.progress.out',
                  'Log.final.out',
                  'SJ.out.tab']

if 'SeparateSAMold' in env['STAR_PARAMS']:
    star_out_files.append('Chimeric.out.sam')

star_targets = [os.path.join(star_out_dir, f) for f in star_out_files]

star = env.Command(star_targets, 
                   star_src, 
                   [chdir_working_cmd + ' && '+\
                    star_cmd + ' && ' +\
                    chdir_parent_cmd])

env.Precious(star[0]) # do not delete old alignment file untill the new one has been computed

## COUNT AND REPORT MAPPED READS
mappings_file = star[0]
mapped_reads_target = os.path.join(star_out_dir, 'STAR_mapped_reads_count.txt')
mapped_reads_cmd    = '''samtools view -F 4 ${SOURCE} '''\
                      '''| cut -f 1 | sort | uniq | wc -l > $TARGET'''
mapped_reads        = env.Command(mapped_reads_target, mappings_file, mapped_reads_cmd)


Clean('.', star_out_dir)

Return('star mapped_reads')
