'''
This script runs the cuffmerge utility from the Cufflinks package.
Software dependencies:
 * Cufflinks v2.2.1

When called from a SConscript it imports the following variables:
 * transcripts_gtf_files
 * cuffmerge_cpus
 * cuffmerge_annotation
 * cuffmerge_genome_fasta

'''

import os, itertools, zipfile

Import('*')

try:
    env = env_cuffmerge.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('GENOME_FASTA', '''The path to the genome single FASTA file, or a directory with one '''
                             '''FASTA file for each chromosome.''', '.')
    vars.Add('TRANSCRIPT_FILES', 'Comma separated list of transcript GTF files', 'sample1_transcripts.gtf,sample2_transcripts.gtf')
    vars.Add('ANNOTATION', 'Genome annotation', 'annotation.gtf')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)
   
    env['TRANSCRIPT_FILES'] = [File(f) for f in env['TRANSCRIPT_FILES'].split(',')]

## COLLECT TRANSCRIPTS.GTF FILES
def WriteAssemblies(target, source, env):
    with open(target[0].path, 'w') as af:
        for gtf in source:
            print>>af, gtf.path

write_af = env.Command('assemblies.txt', 
                       env['TRANSCRIPT_FILES'], 
                       WriteAssemblies)

## MERGE TRANSCRIPTS.GTF FILES
add_options = ' -s ' + env['GENOME_FASTA']
if not env['ANNOTATION'] == '':
    add_options = add_options + ' -g ' + env['ANNOTATION']

cuffmerge_out_dir = 'merged'
cuffmerge_cmd = 'cuffmerge $(-p $CPUS $) -o ${TARGET.dir} ' +\
                add_options + ' ${SOURCE}'
cuffmerge = env.Command(os.path.join(cuffmerge_out_dir, 'merged.gtf'), 
                        [write_af, env['TRANSCRIPT_FILES']], 
                        cuffmerge_cmd)

## AMEND DIRECTORY FROM LOG FILES WHEN CLEANING AND THE DIRECTORY ITSELF
Clean(cuffmerge, os.path.join(cuffmerge_out_dir, 'logs'))
Clean('.', cuffmerge_out_dir)

Return('cuffmerge')
