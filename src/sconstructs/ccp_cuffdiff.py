'''
This SConscript performs differential expression analysis by means of 
Cuffdiff [1]

[1] Trapnell, C. et al. 
    
    Differential analysis of gene regulation at transcript resolution with RNA-seq. 
    
    Nat Biotech 31, 46-53 (2013).

Software dependencies:
 * Cuffdiff v2.2.1

Variables and functions required to export when calling from a SConscript:
 * cuffdiff_cpus
 * cuffdiff_conditions
 * cuffdiff_quantities
 * cuffdiff_annotation

'''

import os, re, collections 

Import('*')

def get_node_sample(node, samples):
    s = None
    for sample in samples:
        if re.match('.*'+sample+'.*', node.path):
            s = sample
            break
    return s

def get_node_condition(node, conditions):
    condition = None
    for c,samples in conditions.iteritems():
        s = get_node_sample(node, samples)
        if s:
            condition = c
            break
    return condition

def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def get_matching_nodes(nodelist, rexpression):
    files = []
    for node in flatten(nodelist):
        if re.match(rexpression, node.path):
            files.append(node)
    return files

def format_cuffdiff_sample_arguments(conditions, cuffdiff_quantities):
    labels = []
    node_lists = ''
    for condition in sorted(conditions.keys()):
        labels.append(condition)
        nodes = []
        for sample in conditions[condition]:
            nodes.extend(get_matching_nodes(cuffdiff_quantities, '.*'+sample+'.*'))
        node_lists = node_lists + ','.join([n.path for n in nodes]) + ' '
    return {'LABELS':','.join(labels), 'SAMS':node_lists}
 
cuffdiff_parameters = {}

try:
    env = env_diffexp.Clone()
    cuffdiff_labels = [get_node_condition(n, env['CONDITIONS']) for n in env['EXPRESSION_FILES']]
    cuffdiff_parameters = format_cuffdiff_sample_arguments(env['CONDITIONS'], env['EXPRESSION_FILES'])
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('QUANTITIES', 'Comma separated list of condition replicates alignment files (SAM/BAM) or Cuffquant files (CXB); @ character separates condition lists', 
             'cond1_1.cxb,cond1_2.cxb@cond2_1.cxb,cond2_2.cxb') 
    vars.Add('ANNOTATION', 'Gene annotation (Ensembl GFF)', 'merged.gtf')
    vars.Add('CONDITIONS', '''Conditions to test differential expression comma separated list ordered as the QUANTITIES parameter''', 
             'CONDITION_1,CONDITION_2')
    vars.Add('EXTRA_PARAMS', 
             'Parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA '\
             ' --multi-read-correct --max-bundle-frags 9999999', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

    CPUS = env['CPUS']
    QUANTITIES = ' '.join(env['QUANTITIES'].split('@'))

    cuffdiff_parameters = {'LABELS': CONDITIONS, 'SAMS': QUANTITIES}

cuffdiff_cmd = 'cuffdiff -q $(--no-update-check -p $CPUS $) -o ' + Dir('.').path + ' ' +\
               '$EXTRA_PARAMS -L ' + cuffdiff_parameters['LABELS'] + ' ${SOURCES[0]} ' +\
               cuffdiff_parameters['SAMS']

cuffdiff_targets = ['isoforms.fpkm_tracking', 'genes.fpkm_tracking', 
                    'cds.fpkm_tracking', 'tss_groups.fpkm_tracking',
                    'isoforms.count_tracking', 'genes.count_tracking', 
                    'cds.count_tracking', 'tss_groups.count_tracking',
                    'isoforms.read_group_tracking', 'genes.read_group_tracking',
                    'cds.read_group_tracking', 'tss_groups.read_group_tracking',
                    'isoform_exp.diff', 'gene_exp.diff',
                    'tss_group_exp.diff', 'cds_exp.diff',
                    'splicing.diff', 'cds.diff',
                    'promoters.diff', 'read_groups.info', 'run.info']
cuffdiff = env.Command(cuffdiff_targets,  
                       [env['ANNOTATION'], env['EXPRESSION_FILES']], 
                       cuffdiff_cmd)

results = {'GENE_EXPRESSION': cuffdiff[9], #genes.read_group_tracking
           'META': cuffdiff[19], #read_groups.info 
           'GENE_DIFF_TESTS': cuffdiff[13], #gene_exp.diff 
           'ALL': cuffdiff}

Return('results')
