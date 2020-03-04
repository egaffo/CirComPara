'''
This SConscript impleemnts the building of a genome index for the STAR read aligner
Import:
    * env_index_star environment including CPUS, GENOME, and EXTRA_PARAMS variables.
                       The GENOME environment variable must report absolute paths of 
                       the FASTA files
Returns:
    * ['chrLength.txt', 'chrNameLength.txt', 
       'chrName.txt', 'chrStart.txt', 
       'Genome', 'genomeParameters.txt', 
       'SA', 'SAindex', 
       'Log.out']
'''

import os

Import('*')

try:
    env = env_index_star
except NameError as ne:
    print 'ccp_index_star.scons: command line execution'
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
            'Comma separated', '')
    vars.Add('EXTRA_PARAMS', 'Extra parameters for STAR', '')
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash', variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)
    
SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src')
    
CPUS        = env['CPUS']
GENOME      = env['GENOME'].split(',')
EXTRA_PARAMS= env['EXTRA_PARAMS']

source  = [File(f).abspath for f in GENOME]
target_basename = '_'.join([os.path.splitext(os.path.basename(f))[0] for f in GENOME])
index_files = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'Genome', 
               'genomeParameters.txt', 'SA', 'SAindex', 'Log.out']
target  = [os.path.join(target_basename, f) for f in index_files]
command = '''cd ${TARGETS[0].dir} && STAR --runMode genomeGenerate $(--runThreadN ''' + CPUS +\
          ''' $) --genomeFastaFiles ''' + ' '.join(source) +\
          ''' --genomeDir . ''' + EXTRA_PARAMS 
index = env.Command(target, source, [command + ''' && cd ''' + Dir('#').abspath])

Clean('.', target_basename)

Return('index')
