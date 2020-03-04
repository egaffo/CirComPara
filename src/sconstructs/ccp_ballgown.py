import os

Import('*')

try:
    env = env_diffexp.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('META', '', '')
    vars.Add('EXPRESSION_FILES', 'A named list of samples ctab files', '')
    vars.Add('EXTRA_PARAMS', '', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    #TODO: 
    #env.Replace(EXPRESSION_FILES = from list to dictionary)

ballgown_targets = 'ballgown.rda'

ballgown_sources = [env['META'], 
                    env['EXPRESSION_FILES'].values()]


ballgown_ctabs =  ",".join(["'" + sample + "'" + "='" + sample_link[0].dir.path + "'" \
                            for sample,sample_link in env['EXPRESSION_FILES'].items()])

ballgown_R_cmd = '''library(data.table);'''\
                 '''meta <- data.frame(unique(fread('${SOURCES[0].path}')[,'''\
                 '''.(id = sample, sample, condition)]), row.names = 'sample');'''\
                 '''library(ballgown);'''\
                 '''samps <- c(''' + ballgown_ctabs + '''); '''\
                 '''bg = ballgown(samples = samps, meas='all', '''\
                 '''              pData = meta[names(samps), , drop = F]);'''\
                 '''save(bg, file='${TARGETS[0]}')'''

ballgown_cmd = '''Rscript -e "''' + ballgown_R_cmd + '''"'''

ballgown = env.Command(ballgown_targets, 
                       ballgown_sources,
                       ballgown_cmd)

results = {'GENE_EXPRESSION': ballgown, 
           'META': env['META'], 
           'GENE_DIFF_TESTS': ballgown, 
           'ALL': ballgown}


Return('results')
