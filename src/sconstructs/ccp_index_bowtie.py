'''
This SConscript impleemnts the building of a genome index for the Bowtie1 read aligner
Import:
    * env_index_bowtie1 environment including CPUS, GENOME, and EXTRA_PARAMS variables.
                       The GENOME environment variable must report absolute paths of 
                       the FASTA files
Returns:
    * ['.1.ebwt', '.2.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']
'''

import os

def SymLink(target, source, env):
    os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

Import('*')

try:
    env = env_index_bowtie.Clone()
except NameError as ne:
    print 'ccp_index_bowtie.scons: command line execution'
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
            'Comma separated', '')
    vars.Add('BOWTIE_EXTRA_PARAMS', 'Extra parameters for bowtie-build', '')
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash', variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)
    
GENOME      = env['GENOME'].split(',')

source  = [File(f).abspath for f in GENOME]
target_basename = '_'.join([os.path.splitext(os.path.basename(f))[0] for f in GENOME])
index_suffixes = ['.1.ebwt', '.2.ebwt', 
                '.3.ebwt', '.4.ebwt',
                '.rev.1.ebwt', '.rev.2.ebwt']
target  = [target_basename + suffix for suffix in index_suffixes]
command = '''bowtie-build -f --seed 1 '''+\
          '''$BOWTIE_EXTRA_PARAMS ''' + ','.join(source) +\
          ''' ${TARGETS[0].dir}''' + os.path.sep + target_basename 
index = env.Command(target, source, command)

link_genome_fasta = env.Command(os.path.join(str(File(index[0]).dir), '${SOURCES[0].file}'), 
                                File(GENOME), 
                                SymLink)

Return('index')
