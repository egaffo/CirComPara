'''
This SConscript runs the cuffquant utility from the Cufflinks software suite.
It quantifies the expression from BAM files.

Software dependencies:
 * cuffquant
When called from a SConscript it imports the following variables:
 * cuffquant_alignments
 * cuffquant_annotation
 * cuffquant_genome_fasta
 * cuffquant_cpus

Returns:
 [al_1.cbx, al_2.cbx]

'''
import os, itertools, zipfile

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_cuffquant.Clone()

except NameError as ne:
    if ne == 'cuffquant_extra_params':
        print 'Warning: ' + ne + ' not defined. Setting no extra paramers as default'
        EXTRA_PARAMS = ''
    else:
        vars = Variables('vars.py')
        vars.Add('CPUS', 'Max parallel jobs to execute', '4')
        vars.Add('GENOME_FASTA', '', '')
        vars.Add('ALIGNMENTS', 'Comma separated list of BAM files', 'al_1.bam,al_2.bam')
        vars.Add('ANNOTATION', 'Annotation GTF file such as from cuffmerge', 'merged.gtf')
        vars.Add('EXTRA_PARAMS', 
                 'Parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA '\
                 ' --multi-read-correct --max-bundle-frags 9999999', '')
        
        env = Environment(ENV=os.environ,variables=vars)
    
        Help(vars.GenerateHelpText(env))
        unknown = vars.UnknownVariables()
        if unknown:
            print "Unknown variables:", unknown.keys()
            Exit(1)
    
        # These are the variables given from the command line when the SConscript is called
        # standalone
        env['ALIGNMENTS'] = [File(f) for f in env['ALIGNMENTS'].split(',')] 

## MERGE TRANSCRIPTS.GTF FILES
add_options = ' --seed 0 ' + env['EXTRA_PARAMS']

cuffquant_cmd = 'cuffquant -q $(--no-update-check -p $CPUS '''\
                '$) -o ${TARGET.dir} ' + add_options +\
                ' ${SOURCES[0]} ${SOURCES[1]}'

samples_abundances = []
for alignment in env['ALIGNMENTS']:
    cuffquant = env.Command(os.path.join('${SOURCES[1].filebase}', 'abundances.cxb'), 
                            [env['ANNOTATION'], alignment], 
                            cuffquant_cmd)
    samples_abundances.append(cuffquant)

Return('samples_abundances')
