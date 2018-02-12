import os 

Import('*')

try:
    env = env_gffcompare.Clone()
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


## GFFCOMPARE (actually, will use CUFFCOMPARE by now since gffcompare
## v0.10.2 and v0.10.3 gives segmentation fault)
if not env['ANNOTATION'] == '':
    
    gffcompare_cmd = "cuffcompare -r ${SOURCES[0]} -R -o "\
                     "${TARGETS[0].dir}/strtcmp -Q ${SOURCES[1]}"
    gffcompare_target = ['strtcmp.combined.gtf', # strtcmp.annotated.gtf, (with gffcompare)
                         'strtcmp.loci', 
                         'strtcmp.stats', 
                         'strtcmp.tracking']
    gffcompare_src = [env['ANNOTATION'], env['MERGED_TRANSCRIPTS']]
    gffcompare = env.Command(gffcompare_target, 
                             gffcompare_src,
                             gffcompare_cmd)

Return('gffcompare')
