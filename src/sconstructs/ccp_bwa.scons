import os

Import('*')

try:
    # these variables can be passed with 'exports' when calling this SConscript
    # from another SConscript
    env             = bwa_env.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('BWA_INDEX', 'The BWA index', '')
    vars.Add('BWA_PARAMS', 'The BWA extra parameters', '')
    vars.Add('READS', 'Input reads. If paired-end, use a comma separated list', 'reads.fa')
    vars.Add('SAMPLE', 'Name of the sample', '')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

#    CPUS            = env['CPUS']
#    BWA_INDEX       = env['BWA_INDEX'] 
    env['READS']     = env['READS'].split(',')
#    SAMPLE          = env['SAMPLE']
#    BWA_PARAMS      = env['BWA_PARAMS']

out_dir = 'bwa_out'

## Single-end reads case
read_bwa_par = ' ${SOURCES[0].abspath} '
## Adjust for paired-end reads
if len(env['READS']) == 2:
    read_bwa_par = read_bwa_par + ' ${SOURCES[1].abspath}'

bwa_target = os.path.join(out_dir, env['SAMPLE'] + '_bwa.sam.gz')

## Map the reads
bwa_cmd = 'bwa mem $(-t $CPUS $) $BWA_PARAMS $BWA_INDEX' +\
          ' ' + read_bwa_par + ' | gzip -c > ${TARGETS[0]}'

## Run the commands
bwa = env.Command(bwa_target,
                  env['READS'], 
                  bwa_cmd)


## COUNT AND REPORT MAPPED READS
mapped_reads_target = os.path.join(out_dir, 'BWA_mapped_reads_count.txt')
mapped_reads_cmd    = '''zcat ${SOURCE} | samtools view -F 4 - '''\
                      '''| cut -f 1 | sort | uniq | wc -l > $TARGET'''
mapped_reads        = env.Command(mapped_reads_target, bwa[0], mapped_reads_cmd)


Return('bwa mapped_reads')
