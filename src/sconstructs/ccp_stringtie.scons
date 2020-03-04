'''
This SConscript runs Stringtie [1] on a single sample.

[1] Pertea, M., Pertea, G.M., Antonescu, C.M., Chang, T.-C., Mendell, J.T., and
    Salzberg, S.L. (2015). StringTie enables improved reconstruction of a
    transcriptome from RNA-seq reads. Nat Biotech 33, 290-295.

Software requirements inherited from the called SConscripts:
 * Stringtie >= v1.3.3b

Variables to include in environment when calling from a SConscript:
    CPUS
    ANNOTATION
    ALIGNMENTS
    SAMPLE
    STRINGTIE_PARAMS
    TOGGLE_TRANSCRIPTOME_RECONSTRUCTION
'''
import os

Import('*')

try:
    env = env_stringtie.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('ANNOTATION', 'The GFF/GTF file with gene annotation (e.g. from Ensembl)', '') 
    vars.Add('ALIGNMENTS', 'The read alignment file in SAM/BAM format', 'sample.bam')
    vars.Add('SAMPLE', 'The sample name', 'sample')
    vars.Add('STRINGTIE_PARAMS', '''Stringtie extra parameters. '''\
             '''F.i. '--rf' assumes a stranded library fr-firststrand, to be '''\
             '''used if dUTPs stranded library were sequenced''', '')
    vars.Add('TOGGLE_TRANSCRIPTOME_RECONSTRUCTION', 'Set True to enable transcriptome '\
	     'reconstruction. Default only quantifies genes and transcripts from the given '\
	     'annotation GTF file', 'False')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

stringtie_targets = [env['SAMPLE'] + '_transcripts.gtf', 
                     env['SAMPLE'] + '_gene_abund.tab']

ballgown_ctabs_dir = 'ballgown_ctabs'

if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
    trx_rec_dir = 'trx_rec'
    ## trx rec if enabled
    if not env['ANNOTATION'] == '':
        env.AppendUnique(STRINGTIE_PARAMS = ['-G', env['ANNOTATION']])

    stringtie_sources = [env['ALIGNMENTS']]
    stringtie_rec_cmd = 'stringtie $( -p $CPUS $) -o ${TARGETS[0]} '\
                    '-l $SAMPLE $STRINGTIE_PARAMS ${SOURCES[0]}'
    stringtie_rec = env.Command(os.path.join(trx_rec_dir, stringtie_targets[0]), 
                            stringtie_sources, 
                            stringtie_rec_cmd)

    ## full merge with original annotation
    if not env['ANNOTATION'] == '':
        
        stringtiemerge_targets = os.path.join(trx_rec_dir,
                                              'stringtiemerge.gtf')
        stringtiemerge_cmd = 'stringtie --merge $( -p $CPUS $) -o ${TARGETS[0]} '\
                             '-G ${SOURCES[1]} ${SOURCES[0]}'
        
        stringtiemerge = env.Command(stringtiemerge_targets,
                                    [stringtie_rec[0], env['ANNOTATION']],
                                    stringtiemerge_cmd)

        env.Replace(ANNOTATION = stringtiemerge[0])
        stringtie_sources.append(env['ANNOTATION'])

    # expression estimate
    stringtie_targets.append(env['SAMPLE'] + '_cov_refs.gtf')
 
    stringtie_targets.extend([os.path.join(ballgown_ctabs_dir, tab) for tab in ['e2t.ctab', 
                                                                                'e_data.ctab',  
                                                                                'i2t.ctab',  
                                                                                'i_data.ctab',  
                                                                                't_data.ctab']]
                            )

    stringtie_cmd = 'stringtie $( -p $CPUS $) -A ${TARGETS[1]} '\
                    '-l $SAMPLE -C ${TARGETS[2]} -b ${TARGETS[3].dir} -e '\
                    '-G $ANNOTATION ${SOURCES[0]} | '\
                    'grep -P "([^\\t]+\\t){6}[+-]\\t" > ${TARGETS[0]}'
                    ## the last piped grep removes from the output the 
                    ## features with unknown strand, which would break 
                    ## TopHat2 runs
    stringtie = env.Command(stringtie_targets,
                            stringtie_sources, 
                            stringtie_cmd)

else:
    if not env['ANNOTATION'] == '':
        if os.path.splitext(env['ANNOTATION'])[1] == '.gz':
            decompressed_annotation = env.Command(os.path.basename(os.path.splitext(env['ANNOTATION'])[0]), 
                                                  env['ANNOTATION'], 
                                                  'zcat $SOURCE > $TARGET')
            env.Replace(ANNOTATION = decompressed_annotation)
    
        env.AppendUnique(STRINGTIE_PARAMS = ['-G', env['ANNOTATION']])
    
        stringtie_targets.append(env['SAMPLE'] + '_cov_refs.gtf')
        env.AppendUnique(STRINGTIE_PARAMS = ['-C', '${TARGETS[2]}'])
     
        stringtie_targets.extend([os.path.join(ballgown_ctabs_dir, tab) for tab in ['e2t.ctab', 
                                                                                    'e_data.ctab',  
                                                                                    'i2t.ctab',  
                                                                                    'i_data.ctab',  
                                                                                    't_data.ctab']]
                                )
        env.AppendUnique(STRINGTIE_PARAMS = ['-b', '${TARGETS[3].dir}'])
    
       
        if not env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
            env.AppendUnique(STRINGTIE_PARAMS = ['-e'])
    
    stringtie_sources = [env['ALIGNMENTS']]
    stringtie_cmd = 'stringtie $( -p $CPUS $) -o ${TARGETS[0]} -A ${TARGETS[1]} '\
                    '-l $SAMPLE $STRINGTIE_PARAMS ${SOURCES[0]}'
    stringtie = env.Command(stringtie_targets, stringtie_sources, stringtie_cmd)
    
    Clean(stringtie, ballgown_ctabs_dir)

return_dict = {'TRANSCRIPTS_GTF': stringtie[0],
               'GENE_EXP': stringtie[1],
               'COV_REFS': '',
               'BALLGOWN_DIR': '',
               'ALL_TARGETS': stringtie}
if not env['ANNOTATION'] == '':
    return_dict['COV_REFS'] = stringtie[2] 
    return_dict['BALLGOWN_DIR'] = File(stringtie[3]).dir

Return('return_dict')

