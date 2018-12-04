'''
COMPUTE BACKSPLICES' LINEAR COUNTS FOR CLR (CIRCULAR TO LINEAR RATIO)
'''

import os, re
Import('*')

def writeListFile(target, source, env):
   
    with open(target[0].path, 'w') as listfile:
        listfile.write('\n'.join([s.abspath for s in source]))

    return None


try:
    env = env_circrna_linexp.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('CIRCRNAS', 'A GTF file with circRNA predictions', '')
    vars.Add('RUNS_DICT', '', '')


    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

env['CIRCOMPARA_HOME'] = env['ENV']['CIRCOMPARA_HOME']


#btmcovstrand = '' ## check if stranded read alignment
#strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")
#if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
#    btmcovstrand = '-s'

## TODO: improve the counting by considering only linear spliced reads
## and linear reads spanning the splice site, while not counting reads
## exactly matching splice site bound (for which we cannot decide 
## whether they belong to the circular or the linear transcript)
#for s in sorted(runs_dict.keys()):
#    ## keep only spliced reads
#    target_1 = "lin_spliced_reads_in_bks.tab"
#    sources_1 = [runs_dict[s]['LINEAR_ALIGNMENTS'], 
#                 circexp['SNP_UNIQUE_CIRC']]
#    cmd_1 = '''samtools view ${SOURCES[0]} | grep "N" | bedtools '''\
#            '''coverage -counts -sorted ''' + btmcovstrand + \
#            ''' -a ${SOURCES[1]} -b stdin > $TARGET '''
#    lin_spliced_reads_in_bks = env.Command(target_1, sources_1, cmd_1)

#    ## keep only linear reads spanning backsplice
#    ## Hint: use bedtools slop to consider bases outside the 
#    ## backsplice start/end sites
#    target_2 = "lin_reads_spanning_bks.tab"
#    sources_2 = sources_1.append(genome_file)
#    cmd_2 = '''samtools view ${SOURCES[0]} | grep -v "N" | '''\
#            '''bedtools coverage -counts -sorted ''' + btmcovstrand + \
#            ''' -a <( bedtools slop -b 1 -i ${SOURCES[1]} ) '''\
#            '''-b stdin > $TARGET '''
#    lin_reads_spanning_bks = env.Command(target_2, sources_2, cmd_2)

#btmcov_sources = [circrna_collect[2], 
#                  [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS'] for 
#                             s in sorted(env['RUNS_DICT'].keys())]]
#btmcov_target = 'bks_linear_counts.tab.gz'
#btmcov_cmd = '''bedtools multicov ''' + btmcovstrand + \
#             ''' -bed ${SOURCES[0]} -bams ${SOURCES[1:]} | gzip -c > ${TARGET}'''
#btmcov = env.Command(btmcov_target, 
#                     btmcov_sources, 
#                     btmcov_cmd)


## use DCC to compute circRNA host gene linear expression
btmcovstrand = '-N' ## check if stranded read alignment
strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")
if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
    btmcovstrand = ''

## convert circRNA GTF to BED as needed by DCC
circ_coords_cmd = '''zcat ${SOURCES[0]} | '''\
                  '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                            '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                            '''([^\\t]+)\\t([^\\t]+)\\tgene_id "([^"]+)";'''\
                            '''/echo -e "\\1\\t\\$$((\\4-1))\\t\\5\\t'''\
                                      '''@\\9@\\t\\6\\t\\7"/e' | '''\
                   '''sed -r 's/@/"/g' >$TARGET '''
circ_coords = env.Command('circrnas.bed', 
                          env['CIRCRNAS'], 
                          circ_coords_cmd)

## make file listing circRNA coordinates for each sample:
## we actually repeat the same circRNA coordinates file
## as it is valid for all samples
env['N_SAMPLES'] = len(env['RUNS_DICT'].keys())
circrna_beds = env.Command('circrna_beds', 
                          [circ_coords for i in range(0, env['N_SAMPLES'])], 
                          writeListFile)

## make file listing linear mapping files for the samples
bams = env.Command('bam_files', 
                   [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS'] for 
                    s in sorted(env['RUNS_DICT'].keys())],
                   writeListFile)

## compute linear counts and set the right names in table header
dcclinexp_sources = [bams, circrna_beds, circ_coords]
dcclinexp_target = 'LinearCount'
dcclinexp_cmd = '''DCC -G ''' + btmcovstrand + \
             ''' -A $GENOME_FASTA -O $TARGET.dir '''\
             '''-B @${SOURCES[0].abspath} -C ${SOURCES[2].abspath} '''\
             '''@${SOURCES[1].abspath} && '''\
             '''sed -i '1s/.*/Chr\\tStart\\tEnd\\t''' +\
             '\\t'.join([s for s in sorted(env['RUNS_DICT'].keys())]) +\
             '''/' $TARGET'''
dcclinexp = env.Command(dcclinexp_target, 
                     dcclinexp_sources, 
                     dcclinexp_cmd)

#for s in sorted(env['RUNS_DICT'].keys()):
#
#    btmcov_sources = [env['CIRCRNAS'], [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS']]
#    btmcov_target = os.path.join(s, 'bks_linear_counts.tab.gz')
#    btmcov_cmd = '''DCC -G ''' + btmcovstrand + \
#                 ''' -A $GENOME_FASTA -B ${SOURCES[1:]} '''\
#                 '''<( zcat ${SOURCES[0]} | cut -f1,4,5 ) | gzip -c > ${TARGET}'''
#    btmcov = env.Command(btmcov_target, 
#                         btmcov_sources, 
#                         btmcov_cmd)
#    

#btmcov_sources = [circrna_collect[2], 
#                  [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS'] for 
#                             s in sorted(env['RUNS_DICT'].keys())]]
#btmcov_target = 'bks_linear_counts.tab.gz'
#btmcov_cmd = '''DCC -G ''' + btmcovstrand + \
#             ''' -A $GENOME_FASTA -B ${SOURCES[1:]} '''\
#             '''<( zcat ${SOURCES[0]} | cut -f1,4,5 ) | gzip -c > ${TARGET}'''
#btmcov = env.Command(btmcov_target, 
#                     btmcov_sources, 
#                     btmcov_cmd)
#
#
Return('dcclinexp')
