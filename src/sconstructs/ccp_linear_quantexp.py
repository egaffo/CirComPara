import os

Import('*')

try:
    env = env_linear_quantexp.Clone()
except NameError:
    print 'ccp_linear_quantexp.scons: cannot import', ne, '. Assume command line script call.'
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('META', '', '')
    vars.Add('LINEAR_ALIGNMENTS', '', '')
    
    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)


## QUANTIFY EXPRESSION
gene_exp = ''
gene_meta = ''

cuffquant_dir = 'cuffquant'
cuffdiff_dir = 'cuffdiff'
htseq_counts_dir = 'htseq_counts'
## the directory where to collect gene expression results
geneexp_dir = 'geneexp'

ballgown_dir = os.path.join(geneexp_dir, 'ballgown')
ballgown_files = []

env['LINEAR_ALIGNMENTS'] = []
for sample,files in env['RUNS_DICT'].iteritems():
    lin_alignments = files['LINEAR_ALIGNMENTS']
    env['LINEAR_ALIGNMENTS'].append(File(lin_alignments))

if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:
    gene_raw_counts_list = []
    trx_raw_counts_list = []
    gene_meta = env['META']
    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
        ## (RE)COMPUTE GENE ABUNDANCES USING MERGED ANNOTATION
        gene_exp = []
        trx_exp  = []
        for sample_name in env['SAMPLES'].keys():
            
            alignment = env['RUNS_DICT'][sample_name]['LINEAR_ALIGNMENTS']

            
            gene_ab_targets = [sample_name + '_transcripts.gtf', 
                               sample_name + '_gene_abund.tab',
                               sample_name + '_cov_refs.gtf']
            
            gene_ab_targets = [os.path.join(geneexp_dir, sample_name, t) for t in gene_ab_targets]

            gene_ab_targets.extend([os.path.join(ballgown_dir, 
                                                 sample_name,
                                                 tab) for tab in ['e2t.ctab', 
                                                                'e_data.ctab',  
                                                                'i2t.ctab',  
                                                                'i_data.ctab',  
                                                                't_data.ctab']]
                                    )
            
            gene_ab_cmd = 'stringtie $( -p $CPUS $) -e -G ${SOURCES[1]} -o ${TARGETS[0]} '\
                          '-A ${TARGETS[1]} -C ${TARGETS[2]} -b ${TARGETS[3].dir} ${SOURCES[0]}' 
            gene_ab = env.Command(gene_ab_targets, 
                                  [alignment, env['ANNOTATION']], 
                                  gene_ab_cmd)
            gene_exp.append(gene_ab[1].abspath)
            trx_exp.append(gene_ab[0].abspath)
            ballgown_files.append(gene_ab[3:])

            ## N.B: this code should be wrapped in a SConscript since it is
            ## used also in ccp_expression.py !
            ## compute read raw counts for StringTie
            SAMPLE = sample_name
            TRANSCRIPTS_GTF = gene_ab[0]
            FASTQC_DATA = env['RUNS_DICT'][SAMPLE]['FASTQC_DATA']
            raw_counts_sources = [TRANSCRIPTS_GTF,
                                  FASTQC_DATA]
            raw_counts_cmd = os.path.join('''get_stringtie_rawcounts.R -g ${SOURCES[0]} '''\
                                     '''-f ${','.join([str(s.abspath) for s in SOURCES[1:]])} '''\
                                     '''-o ${TARGETS[0].dir}''',
                                     SAMPLE + '''_''')
            raw_counts_targets = [os.path.join(SAMPLE + "_" + f) \
                                  for f in ['gene_expression_rawcounts.csv',
                                            'transcript_expression_rawcounts.csv']]
            raw_counts = env.Command(raw_counts_targets, 
                                     raw_counts_sources, 
                                     raw_counts_cmd)
        
            gene_raw_counts_list.append(raw_counts[0])
            trx_raw_counts_list.append(raw_counts[1])

    else:
        ## JUST COLLECT GENE ABUNDANCES ALREADY COMPUTED WITH GIVEN ANNOTATION
        # TODO: use finer way to collect files, by grubbing dictionary results from the runs 
        # gene_exp = []
        # sample_list = []
        # for sample_name in env['SAMPLES'].keys():
        #     gene_exp.append(env['RUNS'][sample_name]['LINEAR_EXPRESSION']['STRINGTIE']['GENE_ABUND'])
        #     sample_list.append(sample_name)
        # sample = ','.join(sample_list)
        # TODO: collect ballgown_files
        gene_exp = get_matching_nodes(env['RUNS'], '.*gene_abund\.tab')
        trx_exp = get_matching_nodes(env['RUNS'], '.*_transcripts\.gtf')
        ballgown_files = get_matching_nodes(env['RUNS'], '.*\.ctab')
        for sample_name in env['SAMPLES'].keys():
            raw_counts = env['RUNS_DICT'][sample_name]['LINEAR_EXPRESSION']['STRINGTIE']['RAW_COUNTS']
            gene_raw_counts_list.append(raw_counts[0])
            trx_raw_counts_list.append(raw_counts[1])

    env.Replace(EXPRESSION_FILES = {'BG_TABS': ballgown_files,
                                    'RAW_COUNTS': gene_raw_counts_list})

   
if ('cufflinks' in env['LINEAR_EXPRESSION_METHODS']) or \
    (env['DIFF_EXP'] and 'cuffdiff' in env['DIFF_EXP']): ## exploit lazy evaluation
    if len(env['SAMPLES'].keys()) > 1:

        env_cuffquant = env.Clone()
        env_cuffquant['EXTRA_PARAMS'] = env['CUFFQUANT_EXTRA_PARAMS']
        env_cuffquant['ALIGNMENTS'] = env['LINEAR_ALIGNMENTS']
        
        cuffquant = SConscript(os.path.join(cuffquant_dir, 'ccp_cuffquant.py'),
                               src_dir = env['SCONSCRIPT_HOME'],
                               variant_dir = cuffquant_dir, duplicate = 0,
                               exports = '''env_cuffquant''')

        Depends(cuffquant, [env['ANNOTATION'], env['LINEAR_ALIGNMENTS']])
        
        ## Normalize expression values: CUFFNORM
        env_cuffnorm = env.Clone()
        env_cuffnorm['QUANT_FILES']= get_matching_nodes(cuffquant, '.*\.cxb')
        env_cuffnorm['EXTRA_PARAMS']= env['CUFFNORM_EXTRA_PARAMS'] #TODO: set strandness when required
        env_cuffnorm['LABELS'] = env['CONDITIONS']
        
        cuffnorm = SConscript(os.path.join(cuffdiff_dir, 'ccp_cuffnorm.py'),
                              src_dir = env['SCONSCRIPT_HOME'],
                              variant_dir = cuffdiff_dir, duplicate = 0,
                              exports = '''env_cuffnorm''')
        
        Depends(cuffnorm, cuffquant)
        
        gene_exp = get_matching_nodes(cuffnorm, ".*genes.read_group_tracking")
        gene_meta = File(get_matching_nodes(cuffnorm, ".*read_groups.info")[0]).abspath
        trx_exp = get_matching_nodes(cuffnorm, ".*isoforms.read_group_tracking")
        gene_raw_counts_list = gene_exp
        trx_raw_counts_list = trx_exp

    else:
        ## case of only one sample
        gene_exp = get_matching_nodes(env['RUNS'], '.*genes.fpkm_tracking')
        gene_meta = env['META']
        #sample = str(env['SAMPLES'].keys()[0])

        trx_exp = get_matching_nodes(env['RUNS'], ".*isoforms.fpkm_tracking")
        gene_raw_counts_list = gene_exp
        trx_raw_counts_list = trx_exp

    env.Replace(EXPRESSION_FILES = gene_exp)

if 'htseq' in env['LINEAR_EXPRESSION_METHODS']:
    ## COMPUTE READ COUNTS FOR EACH ALIGNMENT FILE, ACCORDING TO THE GIVEN ANNOTATION
    htseq_counts = []
    #for algn in get_matching_nodes(env['RUNS'], alignment_matching_regexp):
    #for algn in env['LINEAR_ALIGNMENTS']:
    for sample,files in env['RUNS_DICT'].iteritems():
        env_htseq_count = env.Clone()
        env_htseq_count['SAMPLE'] = sample #os.path.splitext(File(algn).name)[0]
        env_htseq_count['ALIGNMENTS'] = files['LINEAR_ALIGNMENTS'] #File(algn)
        env_htseq_count['STRANDED'] = 'no'

        if '--rna-strandness FR' in env['HISAT2_EXTRA_PARAMS']: 
            env_htseq_count['STRANDED'] = 'yes'

        ## GET READ COUNTS
        htseq_count = SConscript(os.path.join(htseq_counts_dir, 'ccp_htseq_count.py'), 
                                 src_dir = env['SCONSCRIPT_HOME'], 
                                 variant_dir = htseq_counts_dir, duplicate = 0, 
                                 exports = '''env_htseq_count''')

        htseq_counts.append(htseq_count)

    env.Replace(EXPRESSION_FILES = htseq_counts)

    ## TODO: fix gene expression report for htseq_count. 
    gene_meta = env['META']  
    gene_exp = Flatten(htseq_counts)
    trx_exp = gene_exp
    gene_raw_counts_list = gene_exp
    trx_raw_counts_list = gene_exp

## REPORT GENE EXPRESSION
gene_xpr_files_list = env.WriteLinesInTxt(os.path.join(geneexp_dir,
                                          'samples_expression_files.txt'), 
                                  gene_exp)

transcripts_gtf_files_list = env.WriteLinesInTxt(os.path.join(geneexp_dir,
                                          'samples_trxexp_files.txt'), 
                                  trx_exp)

gene_raw_counts_files_list = env.WriteLinesInTxt(os.path.join(geneexp_dir,
                                          'samples_generawc_files.txt'), 
                                  gene_raw_counts_list)

trx_raw_counts_files_list = env.WriteLinesInTxt(os.path.join(geneexp_dir,
                                          'samples_trxrawc_files.txt'), 
                                  trx_raw_counts_list)

gene_exp_analysis_template = os.path.join("$CCP_RMD_DIR", "gene_expression_analysis.Rmd")
gene_exp_analysis_cmd = '''Rscript -e 'results.dir <- dirname("${TARGETS[0].abspath}"); '''\
                        '''meta.file <- "${SOURCES[0].abspath}"; '''\
                        '''gene.xpr.file <- "${SOURCES[1].abspath}"; '''\
                        '''transcripts.gtf.files <- "${SOURCES[2].abspath}"; '''\
                        '''gene_raw_counts_list <- "${SOURCES[3].abspath}"; '''\
                        '''trx_raw_counts_list <-"${SOURCES[4].abspath}"; '''\
                        '''linear.method <- "${ ','.join( LINEAR_EXPRESSION_METHODS ) }"; '''\
                        '''rmarkdown::render(input = "''' + \
                        str(gene_exp_analysis_template) + '''", '''\
                        '''output_file = "$TARGET.abspath", quiet=T,'''\
                        '''intermediates_dir = dirname("${TARGETS[0].abspath}") )' '''

gene_exp_analysis_targets = [os.path.join(geneexp_dir, 
                                          "gene_expression_analysis.html"),
                             os.path.join(geneexp_dir,
                                          "gene_expression_FPKM_table.csv"),
                             os.path.join(geneexp_dir,
                                          "gene_expression_rawcounts_table.csv"),
                             os.path.join(geneexp_dir,
                                         "transcript_expression_rawcounts_table.csv")]

if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:
    gene_exp_analysis_targets.append(os.path.join(geneexp_dir,
                                                  'gene_expression_Nreads_table.csv'))
    gene_exp_analysis_targets.append(os.path.join(geneexp_dir,
                                                  'gene_expression_TPM_table.csv'))

gene_exp_analysis_sources = [gene_meta, 
                             gene_xpr_files_list, 
                             transcripts_gtf_files_list, 
                             gene_raw_counts_files_list,
                             trx_raw_counts_files_list,
                             gene_exp, trx_exp, gene_raw_counts_list, trx_raw_counts_list]
gene_exp_analysis = env.Command(gene_exp_analysis_targets, 
                                gene_exp_analysis_sources, 
                                gene_exp_analysis_cmd)

env['GENE_EXP_ANALYSIS'] = gene_exp_analysis
env['GENE_EXP'] = gene_exp

Clean('.', cuffquant_dir)
Clean('.', cuffdiff_dir)
Clean('.', htseq_counts_dir)
Clean('.', geneexp_dir)
Clean('.', ballgown_dir)

Return('env')
