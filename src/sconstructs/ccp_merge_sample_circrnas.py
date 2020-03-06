'''
'''

import os
Import('*')

try:
    env = env_merge_sample_circrnas.Clone()

except NameError:
    vars = Variables('vars.py')
    vars.Add('SAMPLES', 'List of sample names, comma separated', '')
    vars.Add('RUNS', 'Results of the circRNA methods (comma separated): splicesites.bed for'\
                    'testrealign; circ_candidates.bed for find_circ; ciri.out for CIRI;'\
                    'CIRCexplorer_circ.txt for CIRCexplorer', '')

    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

    env.Replace(SAMPLES, env['SAMPLES'].split(','))
    env.Replace(RUNS, env['RUNS'].split(','))
    SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')
    env['SCONSCRIPT_HOME'] = SCONSCRIPT_HOME

## COLLECT CIRCRNA RESULTS
circRNA_collect_dir = 'merged_samples_circrnas'

samples = sorted(env['SAMPLES']) #samples.keys())
runs = env['RUNS_DICT']
SCONSCRIPT_HOME = env['SCONSCRIPT_HOME']

testrealign_nodes  = []
find_circ_nodes    = []
ciri_nodes         = []
CIRCexplorer_nodes = []
CIRCexplorer2_star_nodes    = []
CIRCexplorer2_bwa_nodes     = []
CIRCexplorer2_segemehl_nodes = []
CIRCexplorer2_tophat_pe_nodes = []
CIRCexplorer2_tophat_nodes = []
dcc_nodes = []
cfinder_nodes = []

## grab circrna methods' results
for sample in samples:
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['TR']
    testrealign_nodes.append(r['CIRCRNAS']) if r else testrealign_nodes.append([''])

    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['FC']
    find_circ_nodes.append(r['CIRCRNAS']) if r else find_circ_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CIRI']
    ciri_nodes.append(r['CIRCRNAS']) if r else ciri_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_ST']
    CIRCexplorer2_star_nodes.append(r['CIRCRNAS']) if r else CIRCexplorer2_star_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_BWA']
    CIRCexplorer2_bwa_nodes.append(r['CIRCRNAS']) if r else CIRCexplorer2_bwa_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_SE']
    CIRCexplorer2_segemehl_nodes.append(r['CIRCRNAS']) if r else CIRCexplorer2_segemehl_nodes.append([''])
    
    ##r = runs[sample]['CIRCULAR_EXPRESSION']['CE2_TH']
    ##CIRCexplorer2_tophat_pe_nodes.append(r['CIRCRNAS']) if r else CIRCexplorer2_tophat_pe_nodes.append([''])
    ##CIRCexplorer2_tophat_pe_nodes.append(get_matching_nodes(runs[sample]['CIRCULAR_EXPRESSION']['CE2_TH'][1])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_TH']
    CIRCexplorer2_tophat_nodes.append(r['CIRCRNAS']) if r else CIRCexplorer2_tophat_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['DCC']
    dcc_nodes.append(r['CIRCRNAS']) if r else dcc_nodes.append([''])
    
    r = runs[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CFINDER']
    cfinder_nodes.append(r['CIRCRNAS']) if r else cfinder_nodes.append([''])

CIRCexplorer_compared = None
ciri_compared = None
find_circ_compared = None
testrealign_compared = None
CIRCexplorer2_star_compared = None
CIRCexplorer2_bwa_compared = None
CIRCexplorer2_segemehl_compared = None
CIRCexplorer2_tophat_pe_compared = None
CIRCexplorer2_tophat_compared = None
dcc_compared = None
cfinder_compared = None

## collect testrealign
if not ''.join([str(f) for f in Flatten(testrealign_nodes)]) == '':

    testrealign_files = ','.join([f.path for f in Flatten(testrealign_nodes)])
    testrealign_labels = ','.join(samples)
    testrealign_comapared_cmd = '''testrealign_compare.R -l ''' + testrealign_labels +\
                             ''' -i ''' + testrealign_files + ''' -o $TARGET'''
    testrealign_comapared_target = os.path.join(circRNA_collect_dir, 'testrealign_compared.csv')
    testrealign_compared = env.Command(testrealign_comapared_target, 
                                    testrealign_nodes, 
                                    testrealign_comapared_cmd)
    Depends(testrealign_compared, testrealign_nodes)

## collect find_circ
if not ''.join([str(f) for f in Flatten(find_circ_nodes)]) == '':

    find_circ_files = ','.join([f.path for f in Flatten(find_circ_nodes)])
    find_circ_labels = ','.join(samples)
    find_circ_comapared_cmd = '''findcirc_compare.R -l ''' + find_circ_labels +\
                             ''' -i ''' + find_circ_files + ''' -o $TARGET'''
    find_circ_comapared_target = os.path.join(circRNA_collect_dir, 'find_circ_compared.csv')
    find_circ_compared = env.Command(find_circ_comapared_target, 
                                    find_circ_nodes, 
                                    find_circ_comapared_cmd)
    Depends(find_circ_compared, find_circ_nodes)
    
## collect CIRI
if not ''.join([str(f) for f in Flatten(ciri_nodes)]) == '':

    ciri_files = ','.join([f.path for f in Flatten(ciri_nodes)])
    ciri_labels = ','.join(samples)
    ciri_comapared_cmd = '''ciri_compare.R -l ''' + ciri_labels +\
                         ''' -i ''' + ciri_files + ''' -o $TARGET'''
    ciri_comapared_target = os.path.join(circRNA_collect_dir, 'ciri_compared.csv')
    ciri_compared = env.Command(ciri_comapared_target, 
                                    ciri_nodes, 
                                    ciri_comapared_cmd)
    Depends(ciri_compared, ciri_nodes)
    
## collect CIRCexplorer2_star
if not ''.join([str(f) for f in Flatten(CIRCexplorer2_star_nodes)]) == '':
    CIRCexplorer2_star_files = ','.join([f.path for f in Flatten(CIRCexplorer2_star_nodes)])
    CIRCexplorer2_star_labels = ','.join(samples)
    CIRCexplorer2_star_comapared_cmd = '''CIRCexplorer_compare.R -l ''' + CIRCexplorer2_star_labels +\
                             ''' -i ''' + CIRCexplorer2_star_files + ''' -o $TARGET'''
    CIRCexplorer2_star_comapared_target = os.path.join(circRNA_collect_dir, 
                                            'CIRCexplorer2_star_compared.csv')
    CIRCexplorer2_star_compared = env.Command(CIRCexplorer2_star_comapared_target, 
                                    CIRCexplorer2_star_nodes, 
                                    CIRCexplorer2_star_comapared_cmd)
    Depends(CIRCexplorer2_star_compared, CIRCexplorer2_star_nodes)

## collect CIRCexplorer2_bwa
if not ''.join([str(f) for f in Flatten(CIRCexplorer2_bwa_nodes)]) == '':
    CIRCexplorer2_bwa_files = ','.join([f.path for f in Flatten(CIRCexplorer2_bwa_nodes)])
    CIRCexplorer2_bwa_labels = ','.join(samples)
    CIRCexplorer2_bwa_comapared_cmd = '''CIRCexplorer_compare.R -l ''' + CIRCexplorer2_bwa_labels +\
                             ''' -i ''' + CIRCexplorer2_bwa_files + ''' -o $TARGET'''
    CIRCexplorer2_bwa_comapared_target = os.path.join(circRNA_collect_dir, 
                                        'CIRCexplorer2_bwa_compared.csv')
    CIRCexplorer2_bwa_compared = env.Command(CIRCexplorer2_bwa_comapared_target, 
                                    CIRCexplorer2_bwa_nodes, 
                                    CIRCexplorer2_bwa_comapared_cmd)
    Depends(CIRCexplorer2_bwa_compared, CIRCexplorer2_bwa_nodes)

## collect CIRCexplorer2_segemehl
if not ''.join([str(f) for f in Flatten(CIRCexplorer2_segemehl_nodes)]) == '':
    CIRCexplorer2_segemehl_files = ','.join([f.path for f in Flatten(CIRCexplorer2_segemehl_nodes)])
    CIRCexplorer2_segemehl_labels = ','.join(samples)
    CIRCexplorer2_segemehl_comapared_cmd = '''CIRCexplorer_compare.R -l ''' + CIRCexplorer2_segemehl_labels +\
                         ''' -i ''' + CIRCexplorer2_segemehl_files + ''' -o $TARGET'''
    CIRCexplorer2_segemehl_comapared_target = os.path.join(circRNA_collect_dir, 
                                        'CIRCexplorer2_segemehl_compared.csv')
    CIRCexplorer2_segemehl_compared = env.Command(CIRCexplorer2_segemehl_comapared_target, 
                                    CIRCexplorer2_segemehl_nodes, 
                                    CIRCexplorer2_segemehl_comapared_cmd)
    Depends(CIRCexplorer2_segemehl_compared, CIRCexplorer2_segemehl_nodes)

## collect CIRCexplorer2_tophat
if not ''.join([str(f) for f in Flatten(CIRCexplorer2_tophat_pe_nodes)]) == '':
    CIRCexplorer2_tophat_nodes = CIRCexplorer2_tophat_pe_nodes

if not ''.join([str(f) for f in Flatten(CIRCexplorer2_tophat_nodes)]) == '':
    CIRCexplorer2_tophat_files = ','.join([f.path for f in Flatten(CIRCexplorer2_tophat_nodes)])    
    CIRCexplorer2_tophat_labels = ','.join(samples)
    CIRCexplorer2_tophat_comapared_cmd = '''CIRCexplorer_compare.R -l ''' + CIRCexplorer2_tophat_labels +\
                         ''' -i ''' + CIRCexplorer2_tophat_files + ''' -o $TARGET'''
    CIRCexplorer2_tophat_comapared_target = os.path.join(circRNA_collect_dir, 
                                        'CIRCexplorer2_tophat_compared.csv')
    CIRCexplorer2_tophat_compared = env.Command(CIRCexplorer2_tophat_comapared_target, 
                                    CIRCexplorer2_tophat_nodes, 
                                    CIRCexplorer2_tophat_comapared_cmd)
    Depends(CIRCexplorer2_tophat_compared, CIRCexplorer2_tophat_nodes)
 
## collect DCC
if not ''.join([str(f) for f in Flatten(dcc_nodes)]) == '':

    dcc_files = ','.join([f.path for f in Flatten(dcc_nodes)])
    dcc_labels = ','.join(samples)
    dcc_comapared_cmd = '''dcc_compare.R -l ''' + dcc_labels +\
                        ''' -i ''' + dcc_files + ''' -o $TARGET'''
    dcc_comapared_target = os.path.join(circRNA_collect_dir, 'dcc_compared.csv')
    dcc_compared = env.Command(dcc_comapared_target, 
                                    dcc_nodes, 
                                    dcc_comapared_cmd)
    Depends(dcc_compared, dcc_nodes)

## collect CircRNA_finder
if not ''.join([str(f) for f in Flatten(cfinder_nodes)]) == '':

    cfinder_files = ','.join([f.path for f in Flatten(cfinder_nodes)])
    cfinder_labels = ','.join(samples)
    cfinder_comapared_cmd = '''cfinder_compare.R -l ''' + cfinder_labels +\
                            ''' -i ''' + cfinder_files + ''' -o $TARGET'''
    cfinder_comapared_target = os.path.join(circRNA_collect_dir, 'cfinder_compared.csv')
    cfinder_compared = env.Command(cfinder_comapared_target, 
                                    cfinder_nodes, 
                                    cfinder_comapared_cmd)
    Depends(cfinder_compared, cfinder_nodes)

## compose a dict with the results
sample_results = {'CIRI': ciri_compared,
                  'FINDCIRC': find_circ_compared,
                  'TESTREALIGN': testrealign_compared,
                  'CIRCEXPLORER2_STAR': CIRCexplorer2_star_compared,
                  'CIRCEXPLORER2_BWA': CIRCexplorer2_bwa_compared,
                  'CIRCEXPLORER2_SEGEMEHL': CIRCexplorer2_segemehl_compared,
                  'CIRCEXPLORER2_TOPHAT': CIRCexplorer2_tophat_compared,
                  'DCC': dcc_compared,
                  'CFINDER': cfinder_compared
                  }

Return('sample_results')
Clean('.', circRNA_collect_dir)
