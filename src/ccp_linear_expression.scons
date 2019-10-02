import os 

Import('*')

try:
    env = env_linear_expression.Clone()
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


## MERGE SAMPLE TRANSCRIPTOMES IF TRANSCRIPTOMES WERE RECONSTRUCTED
cuffmerge_dir = 'new_transcriptome'
transcript_sequences_dir = 'transcript_sequences'
linexp_dir = 'linear_quantexp'

if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
    ## MERGE THE TRANSCRIPTS PREDICTED IN THE VARIOUS SAMPLES
    #transcripts_gtf_files = get_matching_nodes(env['RUNS'], '.*transcripts\.gtf')
    if 'cufflinks' in env['LINEAR_EXPRESSION_METHODS']:
        transcripts_gtf_files = [env['RUNS_DICT'][sample]['LINEAR_EXPRESSION']['CUFFLINKS'][0] for sample in env['RUNS_DICT'].keys()]
    if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:
        transcripts_gtf_files = [env['RUNS_DICT'][sample]['LINEAR_EXPRESSION']['STRINGTIE']['TRANSCRIPTS_GTF'] for sample in env['RUNS_DICT'].keys()]

    if len(transcripts_gtf_files) > 1:        
        if 'cufflinks' in env['LINEAR_EXPRESSION_METHODS']:
            env_cuffmerge = env.Clone()
            env_cuffmerge['TRANSCRIPT_FILES'] = transcripts_gtf_files
            cuffmerge = SConscript(os.path.join(cuffmerge_dir, 'ccp_cuffmerge.scons'),
                                    src_dir = env['SCONSCRIPT_HOME'],
                                    variant_dir = cuffmerge_dir, duplicate = 0,
                                    exports = '''env_cuffmerge''')
            
            env.Replace(ANNOTATION = cuffmerge)

        if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:

            stringtiemerge_list_target = os.path.join(cuffmerge_dir,
                                                      'transcripts_gtf_files_mergelist.txt')
            stringtiemerge_list_cmd = '''echo -ne '${"\\\\n".join([str(s.path) '''\
                                      '''for s in SOURCES])}' > ${TARGETS[0]}''' 
            stringtiemerge_list = env.Command(stringtiemerge_list_target,
                                              transcripts_gtf_files,
                                              stringtiemerge_list_cmd)

            stringtiemerge_targets = [os.path.join(cuffmerge_dir, 'stringtiemerge.gtf')]
            stringtiemerge_cmd = 'stringtie --merge -p $CPUS -o ${TARGETS[0]} '\
            
            if not env['ANNOTATION'] == '':
                stringtiemerge_cmd = stringtiemerge_cmd + ' -G ${SOURCES[1]}' 
            
            stringtiemerge_cmd = stringtiemerge_cmd + ' ${SOURCES[0]}'
            
            stringtiemerge = env.Command(stringtiemerge_targets,
                                        [stringtiemerge_list, env['ANNOTATION']],
                                        stringtiemerge_cmd)
            

            if not env['ANNOTATION'] == '':
                env_gffcompare = env.Clone()
                env_gffcompare['MERGED_TRANSCRIPTS'] = stringtiemerge

                gffcompare = SConscript(os.path.join(cuffmerge_dir, 'ccp_gffcompare.scons'),
                                        src_dir = env['SCONSCRIPT_HOME'],
                                        variant_dir = cuffmerge_dir, duplicate = 0,
                                        exports = '''env_gffcompare''')
                
                stringtiemerge = gffcompare[0]

            env.Replace(ANNOTATION = stringtiemerge)
            
    elif len(transcripts_gtf_files) == 1:
        env.Replace(ANNOTATION = transcripts_gtf_files)


## RETRIEVE TRANSCRIPT SEQUENCES
env_transcript_fastas = env.Clone()
transcript_sequences = SConscript(os.path.join(transcript_sequences_dir, 
                                               'ccp_transcript_fastas.scons'), 
                                  src_dir = env['SCONSCRIPT_HOME'],
                                  variant_dir = transcript_sequences_dir, duplicate = 0,
                                  exports = '''env_transcript_fastas''')

for met in env['LINEAR_EXPRESSION_METHODS']:
    ## QUANTIFY EXPRESSION
    linexp_dir = 'linear_quantexp_' + met
    env_linear_quantexp = env.Clone()
    env_linear_quantexp.Replace(LINEAR_EXPRESSION_METHODS = [met])
    linquantexp = SConscript(os.path.join(linexp_dir, 'ccp_linear_quantexp.scons'),
                                src_dir = env['SCONSCRIPT_HOME'],
                                variant_dir = linexp_dir, duplicate = 0,
                                exports = '''env_linear_quantexp SymLink get_matching_nodes''')
    
    #env['GENE_EXP'] = linquantexp['GENE_EXP']
    #env['GENE_EXP_ANALYSIS'] = linquantexp['GENE_EXP_ANALYSIS']
    #env['EXPRESSION_FILES'] = linquantexp['EXPRESSION_FILES']
    env_linear_quantexp['GENE_EXP'] = linquantexp['GENE_EXP']
    env_linear_quantexp['GENE_EXP_ANALYSIS'] = linquantexp['GENE_EXP_ANALYSIS']
    env_linear_quantexp['EXPRESSION_FILES'] = linquantexp['EXPRESSION_FILES']
    
    ## DIFFERENTIAL LINEAR EXPRESSION
    if env['DIFF_EXP']:
        for de_met in env['DIFF_EXP'].split(','):

            lindiffexp_dir = os.path.join(linexp_dir, 
                                          'linear_diffexp_' + de_met)
            env_linear_diffexp = env_linear_quantexp.Clone()
            env_linear_diffexp['DIFF_EXP'] = de_met
            ## some protocols are fixed:
            ## HTseq-count -> DESeq
            ## Cufflinks -> Cuffdiff
            ## Stringtie -> {Ballgown,DESeq}
            if met == 'htseq':
                if not de_met == 'deseq':
                    continue
            if met == 'cufflinks':
                if not de_met == 'cuffdiff':
                    continue
            if met == 'stringtie':
                if de_met == 'deseq':
                    env_linear_diffexp['EXPRESSION_FILES'] = linquantexp['EXPRESSION_FILES']['RAW_COUNTS']
                elif de_met == 'ballgown':
                    env_linear_diffexp['EXPRESSION_FILES'] = linquantexp['EXPRESSION_FILES']['BG_TABS']
                elif not de_met in ['ballgown', 'deseq']:
                    continue
               
            #env_linear_diffexp['EXPRESSION_FILES'] = linquantexp['EXPRESSION_FILES']
            lindiffexp = SConscript(os.path.join(lindiffexp_dir,
                                                 'ccp_linear_diffexp.scons'),
                                    src_dir = env['SCONSCRIPT_HOME'],
                                    variant_dir = lindiffexp_dir, duplicate = 0,
                                    exports = '''env_linear_diffexp SymLink get_matching_nodes''')
    
            #lindiffexp['DIFFEXP']
            #lindiffexp['REPORT']
            Clean('.', lindiffexp_dir)
    Clean('.', linexp_dir)


Clean('.', cuffmerge_dir)
Clean('.', transcript_sequences_dir)

Return('env')

