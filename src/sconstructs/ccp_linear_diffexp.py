import os, csv, itertools, collections, re, errno
from collections import defaultdict

Import('*')

try:
    env = env_linear_diffexp.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('DIFF_EXP', '', '')
    vars.Add('META', '', '')
    vars.Add('CPUS', '', '')
    vars.Add('EXTRA_PARAMS', '', '')
    vars.Add('SCONSCRIPT_HOME', '', '')
    vars.Add('EXPRESSION_FILES', '', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

conditions  = defaultdict(set)

with open(env['META']) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        conditions[row['condition']].add(row['sample'])

env['SAMPLES'] = Flatten([list(s) for s in conditions.values()])

diffexp = None
report  = None

env['DIFF_EXP'] = env['DIFF_EXP'].lower()

if len(conditions.keys()) > 1:

    for met in env['DIFF_EXP'].split(','):
        env_diffexp = env.Clone()

        ## PREPARE VARIABLES FOR BALLGOWN
        if met == 'ballgown':
            diffexp_sconscript = 'ccp_ballgown.py'
            diffexp_dir = 'ballgown'
            
            ctabs = ['e2t.ctab', 
                     'e_data.ctab', 
                     'i2t.ctab', 
                     'i_data.ctab',
                     't_data.ctab']

            bg_dirs = {}
            for sample in env['SAMPLES']:
                bg_dirs[sample] = []
                for ctab in ctabs:
                
                    rexpr = '.*' + os.path.sep + sample + os.path.sep + '.*' + ctab
                    sample_ctabs = get_matching_nodes(env['EXPRESSION_FILES'], 
                                                      rexpr)
                    sample_link = os.path.join(diffexp_dir, sample, ctab)
                    bg_dirs[sample].extend(env.Command(sample_link,
                                                  sample_ctabs,
                                                  SymLink))

            env_diffexp.Replace(EXPRESSION_FILES = bg_dirs)
            
            report_template = 'ballgown_gene_diffexp.Rmd'

        ## PREPARE VARIABLES FOR CUFFDIFF
        if met == 'cuffdiff':
            diffexp_sconscript = 'ccp_cuffdiff.py'
            diffexp_dir = 'cuffdiff'
            env.Replace(EXTRA_PARAMS = env['CUFFDIFF_EXTRA_PARAMS'])
            env['CONDITIONS'] = conditions
            
            report_template = "cuffdiff_gene_diffexp.Rmd"
            
   
        ## PREPARE VARIABLES FOR DESEQ
        if met == 'deseq':
            diffexp_sconscript = 'ccp_DESeq.py'
            diffexp_dir = 'DESeq2'
            
            report_template = 'DESeq_gene_diffexp.Rmd'
            

        ## COMPUTE DIFFERENTIAL EXPRESSION
        diffexp = SConscript(os.path.join(diffexp_dir, diffexp_sconscript),
                         src_dir = env['SCONSCRIPT_HOME'],
                         variant_dir = diffexp_dir, duplicate = 0,
                         exports = '''env_diffexp''')
        
        if met == 'ballgown':
            Depends(diffexp['GENE_DIFF_TESTS'], bg_dirs.values())
        
        ## MAKE DIFFERENTIAL LINEAR EXPRESSION REPORT
        gene_exp = diffexp['GENE_EXPRESSION']
        gene_meta = diffexp['META']
        gene_diffexp = diffexp['GENE_DIFF_TESTS']
        
        gene_diffexp_analysis_template = os.path.join("$CCP_RMD_DIR", 
        							report_template)
        report_target = [report_template.split('.')[0] + ".html",
                         'DEG_tests_by_contrast.csv']
        report_cmd = '''Rscript -e 'results.dir <- dirname("${TARGETS[1].abspath}"); '''\
        	                '''meta.file <- "${SOURCES[1].abspath}"; '''\
                            '''child.Rmd <- file.path("$CCP_RMD_DIR", "_lindiffexp.Rmd"); '''\
                            '''gene.de.file <- "${SOURCES[0].abspath}"; '''\
                            '''rmarkdown::render(input = "''' + \
                            str(gene_diffexp_analysis_template) + '''", '''\
        	                '''output_file = "${TARGETS[0].abspath}", quiet=T,'''\
        	                '''intermediates_dir = "${TARGETS[0].dir}" )' '''
        report = env.Command(report_target, 
                             [gene_diffexp, gene_meta],
                             report_cmd)

        Clean('.', diffexp_dir)
        
else:
    print "No conditions to contrast in differential linear expression analysis"

results = {'DIFFEXP': diffexp,
           'REPORT': report}

Return('results')
