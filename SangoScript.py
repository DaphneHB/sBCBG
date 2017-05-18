#!/usr/bin/python
# -*- coding: utf-8 -*-

#import shlex
# import subprocess
import os
import time

print "*** /!\ For simulations on Sango, check that interactive=False in LGNeurons.py ***"
print "*** and that the code to be executed has been made executable: chmod +x toto.py ***"

header = '#!/bin/bash \n\n'

execTime = time.localtime()
timeString = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])

print 'Time:', timeString

global i
i = 0 # total number of jobs launched
nbQueuedJobs = 0 # number of jobs currently in the queue

# processes for one parameterization test:
#----------------------------------------- 
def launchOneParameterizedRun(i):
  # checks that there are not too many jobs launched, otherwise waits
  readyToGo = False
  while not readyToGo:
    os.system('squeue -u daphne-heraiz > squeueStatus.txt')

    qf = open('squeueStatus.txt','r')
    nbQueuedJobs = len(qf.readlines())-1
    qf.close()
    
    if nbQueuedJobs < 1000:
      readyToGo = True
      print nbQueuedJobs,"jobs -- Ready to launch job",i
    else:
      print nbQueuedJobs,"jobs -- Wait"
      time.sleep(60)

  # ready to launch a new job:
  IDstring = timeString+'_%05d' %(i)

  print 'Create subdirectory:',IDstring
  os.system('mkdir '+IDstring)
  os.system('cp LGneurons.py '+IDstring+'/')
  os.system('cp testFullBG.py '+IDstring+'/')
  os.system('cp plot_tools.py '+IDstring+'/')
  os.system('cp testChannelBG.py '+IDstring+'/')
  os.system('cp solutions_simple_unique.csv '+IDstring+'/')
  os.system('cp __init__.py '+IDstring+'/')
  os.chdir(IDstring)
  os.system('mkdir log')

  # creation of the modelParams.py file that will correspond to the run at hand
  mltstr = '''#!/apps/free/python/2.7.10/bin/python

# defines the value of the parameters that will be used by testFullbG.py
# generated by sangoScript.py

interactive = False

params = {'nbcpu':    %s,
          'nbCh' :    %d,
          'LG14modelID': %2d,
          'whichTest': %s, # which test was used to generate the log
          'nbMSN':    %f,
          'nbFSI':    %f,
          'nbSTN':    %f,
          'nbGPe':    %f,
          'nbGPi':    %f,
          'nbCSN':    %f,
          'nbPTN':    %f,
          'nbCMPf':   %f,
          'GMSN':     %4.2f,
          'GFSI':     %4.2f,
          'GSTN':     %4.2f,
          'GGPe':     %4.2f,
          'GGPi':     %4.2f, 
          'IeGPe':    %3.1f,
          'IeGPi':    %3.1f,
          'inDegCSNMSN': 100.,
          'inDegPTNMSN':   1.,
          'inDegCMPfMSN':  1.,
          'inDegFSIMSN':  30., # according to Humphries et al. 2010, 30-150 FSIs->MSN
          'inDegMSNMSN':  70., # according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
          'inDegSTNMSN':   0.,
          'inDegGPeMSN':   0.,
          'inDegCSNFSI':  50.,
          'inDegPTNFSI':   1.,
          'inDegSTNFSI':   2.,
          'inDegGPeFSI':  25.,
          'inDegCMPfFSI':  9.,
          'inDegFSIFSI':  15., # according to Humphries et al., 2010, 13-63 FSIs->FSI
          'inDegPTNSTN':  25.,
          'inDegCMPfSTN':  9.,
          'inDegGPeSTN':  25.,
          'inDegCMPfGPe':  9.,
          'inDegSTNGPe':   8.,
          'inDegMSNGPe':2644.,
          'inDegGPeGPe':  25.,
          'inDegMSNGPi':2644.,
          'inDegSTNGPi':   8.,
          'inDegGPeGPi':  23.,
          'inDegCMPfGPi':  9.,
          }
''' %(testedParameters['nbcpu'],testedParameters['nbch'],testedParameters['lg14modelid'],testedParameters['whichTest'],testedParameters['nbmsn'],testedParameters['nbfsi'],testedParameters['nbstn'],testedParameters['nbgpe'],testedParameters['nbgpi'],testedParameters['nbcsn'],testedParameters['nbptn'],testedParameters['nbcmpf'],testedParameters['gmsn'],testedParameters['gfsi'],testedParameters['gstn'],testedParameters['ggpe'],testedParameters['ggpi'],testedParameters['iegpe'],testedParameters['iegpi'])

  print 'Write modelParams.py'
  paramsFile = open('modelParams.py','w')
  paramsFile.writelines(mltstr)
  paramsFile.close()

  # #SBATCH --mem-per-cpu=1G changed for #SBATCH --mem-per-cpu=200M
  slurmOptions = ['#SBATCH --time='+testedParameters['durationH']+':'+testedParameters['durationMin']+':00 \n',
                  '#SBATCH --partition=compute \n',
                  '#SBATCH --mem-per-cpu=1G \n',
                  '#SBATCH --ntasks=1 \n',
                  '#SBATCH --cpus-per-task='+testedParameters['nbcpu']+' \n',
                  '#SBATCH --job-name=sBCBG_'+IDstring+'\n',
                  '#SBATCH --input=none\n',
                  '#SBATCH --output="'+IDstring+'.out" \n',
                  '#SBATCH --error="'+IDstring+'.err" \n',
                  '#SBATCH --mail-user=daphne-heraiz@oist.jp \n',
                  '#SBATCH --mail-type=BEGIN,END,FAIL \n',
                  ]

  moduleUse = ['module use /apps/unit/DoyaU/.modulefiles/ \n']

  moduleLoad = ['module load nest/2.10 \n']

  # write the script file
  print 'Write slurm script file'
  script = open('go.slurm','w')
  script.writelines(header)
  script.writelines(slurmOptions)
  script.writelines(moduleUse)
  script.writelines(moduleLoad)
  #script.writelines('python testFullBG.py \n')
  #script.writelines('python testChannelBG.py \n')
  script.writelines('time srun --mpi=pmi2 python '+testedParameters['whichTest']+'.py \n')
  script.close()

  # execute the script file
  command = 'sbatch go.slurm'
  os.system(command)

  os.chdir('..')

#===============================

# with which additional parameters?
testedParameters={'durationH':    '04',
                  'durationMin':  '00',
                  'nbcpu':        '8',
                  'whichTest':    'testFullBG',
                  'nbch': 1,
                  'lg14modelid':  2,
                  'nbmsn':2644.,
                  'nbfsi':  53.,
                  'nbstn':   8.,
                  'nbgpe':  25.,
                  'nbgpi':  14.,
                  'nbcsn':3000.,
                  'nbptn': 100.,
                  'nbcmpf':  9.,
                  'gmsn':    4.,
                  'gfsi':    1.,
                  'gstn':    1.4,
                  'ggpe':    1.,
                  'ggpi':    1.,
                  'iegpe':  11.,
                  'iegpi':  11.,
                  }

testedParametersIntervals = {}

testedParametersIntervals['lg14modelid']=[2.]
testedParametersIntervals['nbmsn']=[2644.]
testedParametersIntervals['nbfsi']=[  53.]
testedParametersIntervals['nbstn']=[   8.]
testedParametersIntervals['nbgpe']=[  25.]
testedParametersIntervals['nbgpi']=[  14.]
testedParametersIntervals['nbcsn']=[3000.]
testedParametersIntervals['nbptn']=[ 100.]
testedParametersIntervals['nbcmpf']=[  9.]
'''
testedParametersIntervals['gmsn']=[4.]
testedParametersIntervals['gfsi']=[1.]
testedParametersIntervals['gstn']=[1.4]
testedParametersIntervals['ggpe']=[1.]
testedParametersIntervals['ggpi']=[1.]
testedParametersIntervals['iegpe']=[11.]
testedParametersIntervals['iegpi']=[11.]
'''

testedParametersIntervals['gmsn']=[1.,2.,3.,4.,4.5,5.,6.]
testedParametersIntervals['gfsi']=[0.6,0.8,1., 1.2,1.4,1.6,1.8,2.]
testedParametersIntervals['gstn']=[1., 1.1, 1.2,1.3,1.4]
testedParametersIntervals['ggpe']=[0.8,1., 1.1, 1.2,1.4,1.6,1.8]
testedParametersIntervals['ggpi']=[1., 2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,13.]
testedParametersIntervals['iegpe']=[9.,10.,11.,12.,13.,14.,15.]
testedParametersIntervals['iegpi']=[9.,10.,11.,12.,13.,14.,15.]


ftp = open(timeString+'_testedParameter.txt','w')
for k,vlist in testedParametersIntervals.iteritems():
  ftp.writelines(str(k)+' '+str(vlist)+'\n')
ftp.close()

#------------------------------
# recursive exploration of all parameter combinations defined in pdict dictionary
#------------------------------
def recParamExplo(pdict):
  if len(pdict)>0:
    paramK = pdict.keys()[0]
    calldict = pdict.copy()
    del calldict[paramK]
    for v in pdict[paramK]:
      testedParameters[paramK]=v
      recParamExplo(calldict)
  else:
    global i
    launchOneParameterizedRun(i)
    #print i,"->",testedParameters
    i+= 1

recParamExplo(testedParametersIntervals)
