import os, re

Import('*')

try:
    env = env_circular_expression.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('', '', '')

    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

circRNA_collect_dir = 'circRNA_collection'

env_merge_sample_circrnas = env.Clone()
merge_sample_circrnas = SConscript(os.path.join(circRNA_collect_dir, 
					'ccp_merge_sample_circrnas.scons'),
                            src_dir = env['SCONSCRIPT_HOME'], 
                            variant_dir = circRNA_collect_dir, duplicate = 0,
                            exports = '''env_merge_sample_circrnas get_matching_nodes''')


env_circrna_collect = env.Clone()
env_circrna_collect['CSVS'] = merge_sample_circrnas
env_circrna_collect['GTF'] = env['ANNOTATION'] #cuffmerge
circrna_collect = SConscript(os.path.join(circRNA_collect_dir, 
                                          'ccp_collect_circrnas.scons'), 
                            src_dir = env['SCONSCRIPT_HOME'], 
                            variant_dir = circRNA_collect_dir, duplicate = 0,
                            exports = '''env_circrna_collect''')

## COMPUTE BACKSPLICES' LINEAR COUNTS FOR CLR (CIRCULAR TO LINEAR RATIO)
circrna_linexp_dir = 'circrna_linexp'
env_circrna_linexp = env.Clone()
env_circrna_linexp['CIRCRNAS'] = circrna_collect[3] ## unique_circ.gtf.gz

circrna_linexp = SConscript(os.path.join(circrna_linexp_dir, 
                                         'ccp_circrna_linear_expression.py'),
                            src_dir = env['SCONSCRIPT_HOME'],
                            variant_dir = circrna_linexp_dir, duplicate = 0,
                            exports = '''env_circrna_linexp''')
Depends(circrna_linexp, [circrna_collect[2:3], 
                         [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS_IDX'] for s in 
                            env['RUNS_DICT'].keys()]])


## ANALYZE AND REPORT CIRCRNAS 
circrna_analyze_dir = 'circrna_analyze'
env_circrna_analyze = env.Clone()
env_circrna_analyze['META'] = File(env['META']).abspath
env_circrna_analyze['CIRCRNAS'] = circrna_collect[1]
env_circrna_analyze['CIRCGENES'] = circrna_collect[4]
env_circrna_analyze['BKS_LIN_COUNTS'] = circrna_linexp
if int(env_circrna_analyze['MIN_METHODS']) > len(env['CIRCRNA_METHODS']):
    env_circrna_analyze['MIN_METHODS'] = len(env['CIRCRNA_METHODS'])

circrna_analyze = SConscript(os.path.join(circrna_analyze_dir, 
                                          'ccp_analyze_circrnas.scons'),
                            src_dir = env['SCONSCRIPT_HOME'],
                            variant_dir = circrna_analyze_dir, duplicate = 0,
                            exports = '''env_circrna_analyze''')

Clean('.', circRNA_collect_dir)
Clean('.', circrna_analyze_dir)

results = {'SNP_UNIQUE_CIRC': circrna_collect[2]}

Return('results')
