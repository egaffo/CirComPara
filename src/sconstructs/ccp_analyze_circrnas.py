'''
Analyze circRNA predictions and set up a report
'''

import os,re
Import('*')

def write_list_files(target, source, env):
    ''' 
    :param source: a text file with the list of directories
    :param target: must be a 6 element list which will correspond to 
                   - circrnas.gtf.list.txt
                   - circ_to_genes.tsv.list.txt
                   - gene_to_circ.tsv.list.txt
                   - bks_linear_counts.tab.list.txt
                   - meta.csv.list.txt
                   - vars.py.list.txt
    '''
    basedirs = []
    with open(source[0].path, 'r') as dirs:
        for line in dirs:
            basedirs.append(line.strip())

    postfixes = [os.path.join('circular_expression', 'circRNA_collection', 
                              'circrnas.gtf'), 
                 os.path.join('circular_expression', 'circRNA_collection',
                              'combined_circrnas.gtf.gz'), 
                 os.path.join('circular_expression', 'circrna_linexp', 
                              'bks_linear_counts.tab'), 
                 'meta.csv', 
                 'vars.py']
   
    for outfile, postfix in zip(target, postfixes):
       with open(outfile.path, "w") as out:
           for basedir in basedirs:
               out.write(os.path.join(basedir, postfix))
               out.write('\n')

    return None

try:
    env = env_circrna_analyze.Clone()

    ## collect circrnas.gtf files into a filelist
    circrnas_gtf_list_sources = env['CIRCRNAS']
    circrnas_gtf_list_targets = ['circrnas.gtf.list.txt']
    circrnas_gtf_list = env.WriteLinesInTxt(circrnas_gtf_list_targets,
                                            circrnas_gtf_list_sources)
    
    ## collect bks_linear_counts.tab files into a filelist
    bks_linear_counts_list_sources = env['BKS_LIN_COUNTS']
    bks_linear_counts_list_targets = ['bks_linear_counts.tab.list.txt']
    bks_linear_counts_list = env.WriteLinesInTxt(bks_linear_counts_list_targets,
                                                 bks_linear_counts_list_sources)

except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('CIRCRNAS', 
             'A GTF file with circRNA predictions or a text file with a'\
             'list of circrnas.gtf file paths, one per row', 
             'circrnas.gtf.list.txt')
    #vars.Add('CIRC2GENES', 
    #         'circ_to_genes.tsv file as output of gene_annotation.R or '\
    #         'a text file with a list of circ_to_genes.tsv file paths, '\
    #         'one per row', 
    #         'circ_to_genes.tsv.list.txt')
    #vars.Add('GENE2CIRC', 
    #         'gene_to_circ.tsv file as output of gene_annotation.R or '\
    #         'a text file with a list of gene_to_circ.tsv file paths, '\
    #         'one per row', 
    #         'gene_to_circ.tsv.list.txt')
    vars.Add('CIRCGENES', 
             'A 18 fields GTF file with genes intersections with circRNAs '\
             'or ext file with a list of combined_circrnas.gtf.gz file '\
             'paths, one per row', 
             'combined_circrnas.gtf.gz.list.txt')
    vars.Add('BKS_LIN_COUNTS', 
             'The file with linear gene expression for each circRNA or '\
             'a text file with a list of bks_linear_counts.tab file paths, one per row', 
             'bks_linear_counts.tab.list.txt')
    vars.Add('MIN_METHODS', 
             'Number of methods that commmonly detect a circRNA to '\
             ' define the circRNA as reliable', 
             2)
    vars.Add('MIN_READS', 
             'Number of reads to consider a circRNA as expressed', 
             2)
    vars.Add('CIRC_DIFF_EXP', 
             'Enable circRNA differential expression reporting', 
             'False')
    vars.Add('META', 
             'Table specifying sample condition', 
             'meta.csv.list.txt')
    vars.Add('VARS', 
             'vars.py or a text file with a list of vars.py file paths, '\
             'one per row', 
             'vars.py.list.txt')
    vars.Add('CCP_DIRS', 
             'A text file with a list of directories of CirComPara analysis, '\
             'one per row. Usin this option will automatically '\
             'set CIRCRNAS, CIRC2GENES, GENE2CIRC, BKS_LIN_COUNTS, '\
             'VARS, and META parameters by collecting the appropriate files '\
             'form the listed directories. This will ease the merging of '\
             'several runs of CirComPara. N.B: it assumes the meta.csv '\
             'file is named as such and within the CirComPara run directory.', 
             '')

    env = Environment(ENV = os.environ, 
                      SHELL = '/bin/bash',
                      variables = vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)


    #from ccp_utils import writeLines    
    #WriteLinesInTxt = Builder(action = writeLines, suffix = '.txt')
    #env.Append(BUILDERS = {'WriteLinesInTxt' : WriteLinesInTxt})

    env['CIRCOMPARA_HOME'] = env['ENV']['CIRCOMPARA_HOME']

    circrnas_gtf_list = File(env['CIRCRNAS'])
    bks_linear_counts_list = File(env['BKS_LIN_COUNTS'])
    #env['CIRCGENES'] = [File(env['CIRC2GENES']), 
    #                    File(env['GENE2CIRC'])]

    if env['CIRC_DIFF_EXP'].lower() == 'true':
        env['CIRC_DIFF_EXP'] = True
    else:
        env['CIRC_DIFF_EXP'] = False

    WriteListFiles = Builder(action = write_list_files)
    env.Append(BUILDERS = {'WriteListFiles' : WriteListFiles})

    if env['CCP_DIRS'] != '':
        make_lists_targets = ['circrnas.gtf.list.txt', 
                              'combined_circrnas.gtf.gz.list.txt', 
                              'bks_linear_counts.tab.list.txt', 
                              'meta.csv.list.txt', 
                              'vars.py.list.txt']

        ## generate source file path list text files
        make_lists = env.WriteListFiles(make_lists_targets, env['CCP_DIRS'])
        
        circrnas_gtf_list = make_lists[0] #env['CIRCRNAS']
        env['CIRCGENES'] = make_lists[1] #gene annotation
        bks_linear_counts_list = make_lists[2] #env['BKS_LIN_COUNTS']
        env['META'] = make_lists[3]
        env['VARS'] = make_lists[4] 

    ## merge gene annotation 
    merge_gene_annotation_targets = ['circ_to_genes.tsv', 
                                     'gene_to_circ.tsv']
    merge_gene_annotation_cmd = '''gene_annotation.R -c $SOURCE -o ${TARGETS[0].dir}'''
    merge_gene_annotation = env.Command(merge_gene_annotation_targets,
                                   env['CIRCGENES'],
                                   merge_gene_annotation_cmd)
    env['CIRCGENES'] = [merge_gene_annotation[0], #env['CIRC2GENES']
                        merge_gene_annotation[1]] #env['GENE2CIRC']
               

## get circrna raw expression
circrnas_xpr_sources = [circrnas_gtf_list, 
                        bks_linear_counts_list]

circrnas_xpr_targets = ['unfiltered_circrnas.csv', 
                        'ccp_circrna_raw_xpr.csv', 
                        'ccp_circrna_n_methods.csv', 
                        'ccp_circrna_methods.csv',
                        'ccp_bks_linexp.csv'] 
                        ## TODO: add each method expression matrix to targets

circrnas_xpr_command = '''ccp_circrna_expression.R '''\
                       '''-c ${SOURCES[0]} -l ${SOURCES[1]} '''\
                       '''-r ${MIN_READS} -m ${MIN_METHODS} '''\
                       '''-o ${TARGETS[0].dir}'''

circrnas_xpr = env.Command(circrnas_xpr_targets,
                           circrnas_xpr_sources,
                           circrnas_xpr_command)

## make report html
circrnas_analysis_cmd = '''Rscript -e 'results.dir <- dirname("$TARGET.abspath"); '''\
                        '''circrnas.gtf.file <- "${SOURCES[0].abspath}"; '''\
                        '''circ_to_genes.file <- "${SOURCES[1].abspath}"; '''\
                        '''gene_to_circ.file  <- "${SOURCES[2].abspath}"; '''\
                        '''ccp_circrna_raw_xpr.csv.file <- "${SOURCES[3].abspath}"; '''\
                        '''bks_linear_counts.tab.gz.file <- "${SOURCES[4].abspath}"; '''\
                        '''meta_file <- "${SOURCES[5].abspath}"; '''\
                        '''vars.py.filepath <- "${SOURCES[6].abspath}"; '''\
                        '''ccp_circrna_methods.csv.file <- "${SOURCES[7].abspath}"; '''\
                        '''ccp_circrna_n_methods.csv.file <- "${SOURCES[8].abspath}"; '''\
                        '''min_methods <- ${MIN_METHODS}; '''\
                        '''min_reads <- ${MIN_READS}; '''\
                        '''rmarkdown::render(input = "$CCP_RMD_DIR/circRNAs_analysis.Rmd",'''\
                        '''output_file = "$TARGET.abspath", quiet=T,'''\
                        '''intermediates_dir = dirname("$TARGET.abspath") )' '''

circrnas_analysis_targets = ["circRNAs_analysis.html", 
                             "circRNA_expression_per_sample.csv",
                             "methods_shared_circRNA_counts.csv",
                             "cmet_per_circ.csv",
                             "circRNAs_per_gene.csv",
                             "circRNAs_per_gene_per_sample.csv",
                             "circRNAs_per_gene_per_condition.csv",
                             "reliable.circrna.lin.xpr.csv",
                             "reliable.circrna.clr.csv"]

circrnas_analysis_sources = [circrnas_xpr[0], 
                             env['CIRCGENES'][0], 
                             env['CIRCGENES'][1], 
                             circrnas_xpr[1], 
                             circrnas_xpr[4],
                             #env['PROCESSING_READ_STATS'],
                             env['META'],
                             env['VARS'],
                             circrnas_xpr[3],
                             circrnas_xpr[2]]

circrnas_analysis = env.Command(circrnas_analysis_targets, 
                                circrnas_analysis_sources, 
                                circrnas_analysis_cmd)

Clean(circrnas_analysis, 'Figs')

if env['CIRC_DIFF_EXP']:
    circrnas_diffexp_cmd = '''Rscript -e 'results.dir <- dirname("$TARGET.abspath"); '''\
                           '''meta.file <- "${SOURCES[0].abspath}"; '''\
                           '''circrnas.per.sample.file <- "${SOURCES[1].abspath}"; '''\
                           '''rmarkdown::render(input = "$CCP_RMD_DIR/circrnas_diffexp.Rmd",'''\
                           '''output_file = "$TARGET.abspath", quiet=T,'''\
                           '''intermediates_dir = dirname("$TARGET.abspath") )' '''
    
    circrnas_diffexp_targets = "circrnas_diffexp.html"
    circrnas_diffexp_sources = [File(env['META']).abspath, circrnas_analysis[1]]
    circrnas_diffexp = env.Command(circrnas_diffexp_targets, 
                                   circrnas_diffexp_sources, 
                                   circrnas_diffexp_cmd)
      
Return('circrnas_analysis')
