'''
This SConscript build genome indexes for various read aligners.
Currently supported are:
    * HISAT2
    * BWA
    * BOWTIE2
    * SEGEMEHL
    * STAR
'''

import os

Import('*')
try:
    env = env_build_indexes.Clone()
except NameError as ne:
    print 'ccp_build_indexes.scons: command line execution.'

    vars = Variables('vars.py')
    vars.Add('INDEXES', 'A comma separated list of read aligner programs for which'\
                        ' the genome index has to be built', 
             'HISAT2,BWA,BOWTIE2,SEGEMEHL,STAR')
    vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
                       'Comma separated.', '')
    vars.Add('HISAT2_EXTRA_PARAMS', 'Extra parameters for htseq2-build', '')
    vars.Add('BWA_EXTRA_PARAMS', 'Extra parameters for bwa index', '')
    vars.Add('BOWTIE2_EXTRA_PARAMS', 'Extra parameters for bowtie2-build', '')
    vars.Add('SEGEMEHL_EXTRA_PARAMS', 'Extra parameters for segemehl.x -x', '')
    vars.Add('STAR_EXTRA_PARAMS', 'Extra parameters for STAR --genomeGenerate', '')
    vars.Add('BOWTIE_EXTRA_PARAMS', 'Extra parameters for bowtie-build', '')
    
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)


SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

## SET DEFAULT VALUES FOR (OPTIONAL) EXTRA PARAMETERS
env.SetDefault(HISAT2_EXTRA_PARAMS  = '')
env.SetDefault(BWA_EXTRA_PARAMS     = '')
env.SetDefault(BOWTIE2_EXTRA_PARAMS = '')
env.SetDefault(SEGEMEHL_EXTRA_PARAMS = '')
env.SetDefault(STAR_EXTRA_PARAMS    = '')
env.SetDefault(BOWTIE_EXTRA_PARAMS = '')

indexes_dir = 'indexes'
indexes = {}

if 'HISAT2' in env['INDEXES'].split(','):
    
    env_index_hisat2 = env
    env_index_hisat2['CPUS'] = env['CPUS']
    env_index_hisat2['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_hisat2['EXTRA_PARAMS'] = env['HISAT2_EXTRA_PARAMS']

    index_dir = os.path.join(indexes_dir, 'hisat2')
    index = SConscript(os.path.join(index_dir, 'ccp_index_hisat2.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_hisat2 ''')

    indexes['HISAT2'] = index
    Clean('.', index_dir)

if 'SEGEMEHL' in env['INDEXES'].split(','):
    
    env_index_segemehl = env
    env_index_segemehl['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_segemehl['EXTRA_PARAMS'] = env['SEGEMEHL_EXTRA_PARAMS']

    index_dir = os.path.join(indexes_dir, 'segemehl')
    index = SConscript(os.path.join(index_dir, 'ccp_index_segemehl.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_segemehl ''')

    indexes['SEGEMEHL'] = index
    Clean('.', index_dir)

if 'BWA' in env['INDEXES'].split(','):
    
    env_index_bwa = env
    env_index_bwa['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_bwa['EXTRA_PARAMS'] = env['BWA_EXTRA_PARAMS']
    #env_index_hisat2['INDEXING_ALGORITHM'] = env['INDEXING_ALGORITHM']

    index_dir = os.path.join(indexes_dir, 'bwa')
    index = SConscript(os.path.join(index_dir, 'ccp_index_bwa.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_bwa ''')
    indexes['BWA'] = index
    Clean('.', index_dir)


if 'BOWTIE2' in env['INDEXES'].split(','):
    
    env_index_bowtie2 = env
    env_index_bowtie2['CPUS'] = env['CPUS']
    env_index_bowtie2['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_bowtie2['EXTRA_PARAMS'] = env['BOWTIE2_EXTRA_PARAMS']

    index_dir = os.path.join(indexes_dir, 'bowtie2')
    index = SConscript(os.path.join(index_dir, 'ccp_index_bowtie2.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_bowtie2 ''')
    indexes['BOWTIE2'] = index
    Clean('.', index_dir)


if 'BOWTIE' in env['INDEXES'].split(','):
    
    env_index_bowtie = env.Clone()
    env_index_bowtie['CPUS'] = env['CPUS']
    env_index_bowtie['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_bowtie['EXTRA_PARAMS'] = env['BOWTIE_EXTRA_PARAMS']

    index_dir = os.path.join(indexes_dir, 'bowtie')
    index = SConscript(os.path.join(index_dir, 'ccp_index_bowtie.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_bowtie ''')
    indexes['BOWTIE'] = index
    Clean('.', index_dir)


if 'STAR' in env['INDEXES'].split(','):
    
    env_index_star = env
    env_index_star['CPUS'] = env['CPUS']
    env_index_star['GENOME'] = ','.join([File(f).abspath for f in env['GENOME'].split(',')])
    env_index_star['EXTRA_PARAMS'] = env['STAR_EXTRA_PARAMS']

    index_dir = os.path.join(indexes_dir, 'star')
    index = SConscript(os.path.join(index_dir, 'ccp_index_star.py'), 
                            src_dir = SCONSCRIPT_HOME, 
                            variant_dir = index_dir, duplicate = 0,
                            exports = '''env_index_star ''')
    indexes['STAR'] = index
    Clean('.', index_dir)

Clean('.', indexes_dir)

Return('indexes')
