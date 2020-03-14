'''
This SConscript impleemnts the building of a genome index for the BWA read aligner.

Import:
    * env_index_bwa environment including INDEXING_ALGORITHM, GENOME, and EXTRA_PARAMS variables.
                       The GENOME environment variable must report absolute paths of 
                       the FASTA files
Returns:
    * ['.amb', '.ann', '.bwt', '.pac', '.sa']
'''

import os

Import('*')

try:
    env = env_index_bwa
except NameError as ne:
    print 'ccp_index_bwa.scons: command line execution'
    vars = Variables('vars.py')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
            'Comma separated', '')
    vars.Add('EXTRA_PARAMS', 'Extra parameters for segemehl', '')
    vars.Add('INDEXING_ALGORITHM', 'The algorithm to use for index computing', 'bwtsw')
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash', variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)
    
SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

GENOME      = env['GENOME'].split(',')
env.SetDefault(INDEXING_ALGORITHM = 'bwtsw') # set a default value in case it was not defined
INDEXING_ALGORITHM = env['INDEXING_ALGORITHM']
EXTRA_PARAMS= env['EXTRA_PARAMS']

source  = [File(f).abspath for f in GENOME]
target_basename = '_'.join([os.path.splitext(os.path.basename(f))[0] for f in GENOME])
index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
target  = [target_basename + suffix for suffix in index_suffixes]
command = '''bwa index -a ''' + INDEXING_ALGORITHM +\
          ''' -p ${TARGETS[0].dir}''' + os.path.sep + '''${TARGETS[0].filebase} ''' +\
          EXTRA_PARAMS + ' ' + ' '.join(source)
index = env.Command(target, source, command)

Return('index')
