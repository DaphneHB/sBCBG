 #!/apps/free/python/2.7.10/bin/python
# -*- coding: utf-8 -*-    
from LGneurons import *
from modelParams import *
from gurneyParams import *

import nest.raster_plot
#import time
import sys
import os
# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population

#------------------------------------------
# Creates the populations of neurons necessary to simulate a BG circuit
#------------------------------------------
def createBG_MC():
  #==========================
  # Creation of neurons
  #-------------------------
  print '\nCreating neurons\n================'

  nbSim['MSN'] = params['nbMSN']
  createMC('MSN',params['nbCh'])

  nbSim['FSI'] = params['nbFSI']
  createMC('FSI',params['nbCh'])

  nbSim['STN'] = params['nbSTN']
  createMC('STN',params['nbCh'])

  nbSim['GPe'] = params['nbGPe']
  createMC('GPe',params['nbCh'])
  for i in range(len(Pop['GPe'])):
    nest.SetStatus(Pop['GPe'][i],{"I_e":params['IeGPe']})

  nbSim['GPi'] = params['nbGPi']
  createMC('GPi',params['nbCh'])
  for i in range(len(Pop['GPi'])):
    nest.SetStatus(Pop['GPi'][i],{"I_e":params['IeGPi']})

  parrot = True # switch to False at your risks & perils...                                                                                                                   
  nbSim['CSN'] = params['nbCSN']
  createMC('CSN',params['nbCh'], fake=True, parrot=parrot)

  nbSim['PTN'] = params['nbPTN']
  createMC('PTN',params['nbCh'], fake=True, parrot=parrot)

  nbSim['CMPf'] = params['nbCMPf']
  createMC('CMPf',params['nbCh'], fake=True, parrot=False)

  print "Number of simulated neurons:", nbSim

#------------------------------------------
# Connects the populations of a previously created multi-channel BG circuit 
#------------------------------------------
def connectBG_MC(antagInjectionSite,antag):
  print "params = ",params
  exit()
  G = {'MSN': params['GMSN'],
       'FSI': params['GFSI'],
       'STN': params['GSTN'],
       'GPe': params['GGPe'],
       'GPi': params['GGPi'],
      }

  print "Gains on LG14 syn. strength:", G

  #-------------------------
  # connection of populations
  #-------------------------
  print '\nConnecting neurons\n================'
  print "**",antag,"antagonist injection in",antagInjectionSite,"**"
  print '* MSN Inputs'
  connectMC('ex','CSN','MSN', 'focused', inDegree= min(params['inDegCSNMSN'],nbSim['CSN']),gain=G['MSN'])
  connectMC('ex','PTN','MSN', 'focused', inDegree= min(params['inDegPTNMSN'],nbSim['PTN']),gain=G['MSN'])
  connectMC('ex','CMPf','MSN','focused', inDegree= min(params['inDegCMPfMSN'],nbSim['CMPf']),gain=G['MSN'])
  connectMC('in','MSN','MSN', 'focused', inDegree= min(params['inDegMSNMSN'],nbSim['MSN']),gain=G['MSN'])
  connectMC('in','FSI','MSN', 'diffuse', inDegree= min(params['inDegFSIMSN'],nbSim['FSI']),gain=G['MSN']) # diffuse ? focused ?                                               
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    print "PB EST LA !!!!!"
    print "params in degree ",params['inDegSTNMSN'],"\nnb simu ",nbSim['STN']
    connectMC('ex','STN','MSN', 'diffuse', inDegree= min(params['inDegSTNMSN'],nbSim['STN']),gain=G['MSN'])
    exit()
  if alpha['GPe->MSN'] != 0:
    print "alpha['GPe->MSN']",alpha['GPe->MSN']
    connectMC('in','GPe','MSN', 'diffuse', inDegree= min(params['inDegGPeMSN'],nbSim['GPe']),gain=G['MSN']) # diffuse ? focused ?                                             

  print '* FSI Inputs'
  connectMC('ex','CSN','FSI', 'focused', inDegree= min(params['inDegCSNFSI'],nbSim['CSN']),gain=G['FSI'])
  connectMC('ex','PTN','FSI', 'focused', inDegree= min(params['inDegPTNFSI'],nbSim['PTN']),gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
    connectMC('ex','STN','FSI', 'diffuse', inDegree= min(params['inDegSTNFSI'],nbSim['STN']),gain=G['FSI'])
  connectMC('in','GPe','FSI', 'focused', inDegree= min(params['inDegGPeFSI'],nbSim['GPe']),gain=G['FSI'])
  connectMC('ex','CMPf','FSI','focused', inDegree= min(params['inDegCMPfFSI'],nbSim['CMPf']),gain=G['FSI'])
  connectMC('in','FSI','FSI', 'diffuse',inDegree= min(params['inDegFSIFSI'],nbSim['FSI']),gain=G['FSI'])

  print '* STN Inputs'
  connectMC('ex','PTN','STN', 'focused',inDegree= min(params['inDegPTNSTN'],nbSim['PTN']), gain=G['STN'])
  connectMC('ex','CMPf','STN','focused',inDegree= min(params['inDegCMPfSTN'],nbSim['CMPf']), gain=G['STN'])
  connectMC('in','GPe','STN', 'focused',inDegree= min(params['inDegGPeSTN'],nbSim['GPe']), gain=G['STN']) # or diffuse, to be in line with the 2008 model?                    

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connectMC('NMDA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('NMDA','STN','GPe','diffuse', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connectMC('in','MSN','GPe','focused', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connectMC('in','GPe','GPe','diffuse', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe']) # diffuse or focused?                                           
    elif antag == 'NMDA':
      connectMC('AMPA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('AMPA','STN','GPe','diffuse', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connectMC('in','MSN','GPe','focused', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connectMC('in','GPe','GPe','diffuse', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connectMC('NMDA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('NMDA','STN','GPe','diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    elif antag == 'GABAA':
      connectMC('ex','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('ex','STN','GPe', 'diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connectMC('ex','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
    connectMC('ex','STN','GPe', 'diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    connectMC('in','MSN','GPe', 'focused',inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
    connectMC('in','GPe','GPe', 'diffuse',inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'All':
      pass
    elif antag == 'NMDA':
      connectMC('in','MSN','GPi',   'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('AMPA','STN','GPi', 'diffuse', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('in','GPe','GPi',   'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connectMC('AMPA','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connectMC('in','MSN','GPi', 'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('in','GPe','GPi', 'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    elif antag == 'AMPA':
      connectMC('in','MSN','GPi',   'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('NMDA','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('in','GPe','GPi',   'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connectMC('NMDA','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'GABAA':
      connectMC('ex','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('ex','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connectMC('in','MSN','GPi', 'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
    connectMC('ex','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
    connectMC('in','GPe','GPi', 'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    connectMC('ex','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])

#------------------------------------------
# Checks that the BG model parameterization defined by the "params" dictionary can respect the electrophysiological constaints (firing rate at rest).
# If testing for a given antagonist injection experiment, specifiy the injection site in antagInjectionSite, and the type of antagonists used in antag.
# Returns [score obtained, maximal score]
# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population
#------------------------------------------

def checkAvgFR_MC(out=dict(),NbTrials=1,ctx_activity=None,showRasters=False,params={},antagInjectionSite='none',antag='',logFileName='',CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1.,shutOffset=False,offTime=1000,simuTime=5000):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = float(offTime)
  simDuration = float(simuTime) # ms
  # nest.SetKernelStatus({"overwrite_files":True}) # Thanks to use of timestamps, file names should now 
                                                   # be different as long as they are not created during the same second

  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  loadLG14params(params['LG14modelID'])

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nprm in necessaryParams:
    if nprm not in params:
      print "Missing parameter:",np 
      exit()

  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()
  connectBG_MC(antagInjectionSite,antag)
  # lets draw the connectivity matrix
  if not ctx_activity is None :
    print "Multi-trial test"
    # ctx input activities
     #-------------------------
    # prepare the firing rates of the inputs for the 5 steps of the experiment
    #-------------------------  
    
    gCSN = CSNFR[1]-CSNFR[0]
    gPTN = PTNFR[1]-PTNFR[0]
    
    actLevels = np.array([ctx_activity])
    NbTrials = len(ctx_activity)
    CSNrate= gCSN * actLevels + np.ones((NbTrials)) * CSNFR[0]
    PTNrate= gPTN * actLevels + np.ones((NbTrials)) * PTNFR[0]
    
    #-------------------------
    # and prepare the lists of neurons that will be affected by these activity changes
    #-------------------------
    ActPop = {'CSN':[()],'PTN':[()]}
    if 'Fake' in globals():
      if 'CSN' in Fake:
        if PActiveCSN==1.:
         ActPop['CSN']=Fake['CSN']
        else:
          for i in chanRange:
            ActPop['CSN'][0] = tuple(rnd.choice(a=np.array(Fake['CSN'][0]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
      else:
        if PActiveCSN==1.:
         ActPop['CSN']=Pop['CSN']
        else:
          for i in chanRange:
            ActPop['CSN'][0] = tuple(rnd.choice(a=np.array(Pop['CSN'][0]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
      if 'PTN' in Fake:
        if PActivePTN==1.:
          ActPop['PTN']=Fake['PTN']
        else:
          for i in chanRange:
            ActPop['PTN'][0] = tuple(rnd.choice(a=np.array(Fake['PTN'][0]),size=int(nbSim['PTN']*PActivePTN),replace=False))
      else:
        if PActivePTN==1.:
          ActPop['PTN']=Pop['PTN']
        else:
          for i in chanRange :
            ActPop['PTN'][0] = tuple(rnd.choice(a=np.array(Pop['PTN'][0]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in chanRange:
          ActPop['CSN'][0] = tuple(rnd.choice(a=np.array(Pop['CSN'][0]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in chanRange :
          ActPop['PTN'][0] = tuple(rnd.choice(a=np.array(Pop['PTN'][0]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else :
    NbTrials = 1
  
  score = 0

  text=[]
  frstr = "#" + str(params['LG14modelID'])+ " , " + antagInjectionSite + ','
  for trial in range(NbTrials) :
    #-------------------------
    # measures
    #-------------------------
    spkDetect={} # spike detectors used to record the experiment
    expeRate={}
  
    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'
  
    for N in NUCLEI:
      # 1000ms offset period for network stabilization
      if shutOffset and timeStep!=0:
        start = trial*simDuration 
        stop = (trial+1)*simDuration
      else :
        start = offsetDuration + trial*(offsetDuration+simDuration) 
        stop = (trial+1)*(offsetDuration+simDuration)

      spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(trial)+'_'+antagStr+N, "to_file": False, 'start': start ,'stop':stop})

      for i in range(len(Pop[N])):
        nest.Connect(Pop[N][i], spkDetect[N])
  
    #-------------------------
    # Simulation
    #-------------------------
    if not ctx_activity is None :
      print '====== Step',trial,'======'
      print 'Channel 0:',CSNrate[0,trial],PTNrate[0,trial]
      
      nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,trial]})
      nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,trial]})
    if shutOffset and timeStep!=0 :
      nest.Simulate(simDuration)
    else :
      nest.Simulate(simDuration+offsetDuration)
  
    s = '----- RESULTS -----'
    print s
    text.append(s+'\n')
    if antagInjectionSite == 'none':
      validationStr = "\n#" + str(params['LG14modelID']) + " , "
      frstr += "none , "
      for N in NUCLEI:
        strTestPassed = 'NO!'
        expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
        if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
          # if the measured rate is within acceptable values
          strTestPassed = 'OK'
          score += 1
          validationStr += N + "=OK , "
        else : 
        # out of the ranges
          if expeRate[N] > FRRNormal[N][1] :
            difference = expeRate[N] - FRRNormal[N][1]
            validationStr += N + "=+%.2f , " % difference
          else :
            difference = expeRate[N] - FRRNormal[N][0]
            validationStr += N + "=%.2f , " % difference
        frstr += '%f , ' %(expeRate[N])
        s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRNormal[N][0])+' , '+str(FRRNormal[N][1])+')'
        print s
        text.append(s+'\n')
    else:
      validationStr = ""
      frstr += str(antag) + " , "
      for N in NUCLEI:
        expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
        if N == antagInjectionSite:
          strTestPassed = 'NO!'
          if expeRate[N] <= FRRAnt[N][antag][1] and expeRate[N] >= FRRAnt[N][antag][0]:
            # if the measured rate is within acceptable values
            strTestPassed = 'OK'
            score += 1
            validationStr += N + "_" + antag + "=OK , "
          else : 
          # out of the ranges
            if expeRate[N] > FRRNormal[N][1] :
              difference = expeRate[N] - FRRNormal[N][1]
              validationStr += N + "_" + antag + "=+%.2f , " % difference
            else :
              difference = expeRate[N] - FRRNormal[N][0]
              validationStr += N + "_" + antag + "=%.2f , " % difference
          
          s = '* '+N+' with '+antag+' antagonist(s): '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRAnt[N][antag][0])+' , '+str(FRRAnt[N][antag][1])+')'
          print s
          text.append(s+'\n')
        else:
          s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz'
          print s
          text.append(s+'\n')
        frstr += '%f , ' %(expeRate[N])
  
    out[trial] = expeRate['GPi']
  
    s = '-------------------'
    print s
    text.append(s+'\n')
    os.system("rm -f log/*.gdf")


  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','a')
  firingRatesFile.writelines(frstr)
  firingRatesFile.close()

  #print "************************************** file writing",text
  #res = open(dataPath+'OutSummary_'+logFileName+'.txt','a')
  res = open(dataPath+'OutSummary.txt','a')
  res.writelines(text)
  res.close()

  validationFile = open("validationArray.csv",'a')
  validationFile.write(validationStr)
  validationFile.close()
  
  #-------------------------
  # Displays
  #-------------------------
  if showRasters and interactive:
    displayStr = ' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
    for N in NUCLEI:
      nest.raster_plot.from_device(spkDetect[N],hist=True,title=N+displayStr)

    nest.raster_plot.show()

  return (score, 5 if antagInjectionSite == 'none' else 1),out

#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGurneyTest(showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag=''):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = 200.
  simDuration = 800. # ms
  loadLG14params(params['LG14modelID'])

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nepa in necessaryParams:
    if nepa not in params:
      print "Missing parameter:",nepa
      exit()

  nbRecord = 2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels

  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()
  connectBG_MC(antagInjectionSite,antag)

  #-------------------------
  # prepare the firing rates of the inputs for the 5 steps of the experiment
  #-------------------------  
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  activityLevels = np.array([[0,0.4,0.4,0.6,0.4], [0.,0.,0.6,0.6,0.6]])  
  
  CSNrate= gCSN * activityLevels + np.ones((5)) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((5)) * PTNFR[0]

  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop = {'CSN':[(),()],'PTN':[(),()]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in range(2):
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in range(2):
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))

  #-------------------------
  # log-related variables
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,5))

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(nbRecord):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'step , '
  for i in range(nbRecord):
    for N in NUCLEI:
      frstr += N+' ('+str(i)+') , '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)

  #----------------------------------
  # Loop over the 5 steps of the test
  #----------------------------------
  for timeStep in range(5):
    #-------------------------
    # measures                                                                                                                                                  
    #-------------------------
    spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels

    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'

    for i in range(nbRecord):
      for N in NUCLEI:
        spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(timeStep)+'_'+antagStr+N, "to_file": True, 'start':offsetDuration + timeStep*(offsetDuration+simDuration),'stop':(timeStep+1)*(offsetDuration+simDuration)})
        nest.Connect(Pop[N][i], spkDetect[i][N])

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    '''
    if 'Fake' in globals():
      if 'CSN' in Fake:
        nest.SetStatus(Fake['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Fake['CSN'][1],{'rate':CSNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      if 'PTN' in Fake:
        nest.SetStatus(Fake['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Fake['PTN'][1],{'rate':PTNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    else:
      nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
      nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
      nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    '''

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    nest.Simulate(simDuration+offsetDuration)

    for i in range(nbRecord):
      print '------ Channel',i,'-------'
      for N in NUCLEI:
        #strTestPassed = 'NO!'
        expeRate[N][i,timeStep] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
        print 't('+str(timeStep)+')',N,':',expeRate[N][i,timeStep],'Hz'
        frstr += '%f , ' %(expeRate[N][i,timeStep])

    strTestPassed = 'YES!'
    if timeStep == 0:
      for i in range(params['nbCh']):
        if expeRate['GPi'][0,timeStep]<FRRNormal['GPi'][0]:
          strTestPassed = 'NO!'
      meanRestGPi = expeRate['GPi'][:,timeStep].mean()
    elif timeStep == 1:
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9:
        strTestPassed = 'NO!'
    elif timeStep == 2:
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][0,timeStep] < expeRate['GPi'][1,timeStep]:
        strTestPassed = 'NO!'
    elif timeStep == 3:
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
    elif timeStep == 4:
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][0,timeStep] < expeRate['GPi'][1,timeStep]:
        strTestPassed = 'NO!'

    if strTestPassed == 'YES!':
      score +=1

    print '------ Result ------'
    print expeRate['GPi'][0,timeStep],'Hz',expeRate['GPi'][1,timeStep],'Hz',strTestPassed
  
    # write measured firing rates in csv file
    frstr+='\n'
    firingRatesFile.writelines(frstr)

    #-------------------------
    # Displays
    #-------------------------
    '''
    if showRasters and interactive:
      for i in range(nbRecord):
        displayStr = ' Channel '+str(i)
        displayStr+=' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
        #for N in NUCLEI:
        for N in ['MSN','STN']:
          #nest.raster_plot.from_device(spkDetect[i][N],hist=True,title=N+displayStr)
          nest.raster_plot.from_device(spkDetect[i][N],hist=False,title=N+displayStr)

      nest.raster_plot.show()
    '''

  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()

  return score,5
  
#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGurneyTestGeneric(trials_dico,ratio=1.5,offsetTime=200,simuTime=800,shuffled=True,xytab=np.arange(0.,1.,0.1),showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag=''):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = float(offsetTime)
  simDuration = float(simuTime) # ms

  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  print params['inDegSTNMSN']
  loadLG14params(params['LG14modelID'])

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nepa in necessaryParams:
    if nepa not in params:
      print "Missing parameter:",nepa
      exit()

  nbRecord = 2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels

  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  print params['inDegSTNMSN']
  createBG_MC()
  print params
  print "entering connection"
  connectBG_MC(antagInjectionSite,antag)
  print params['inDegSTNMSN']
  exit()
  #-------------------------
  # prepare the firing rates of the inputs for the 5 steps of the experiment
  #-------------------------  
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  activityLevels1 = list()
  activityLevels2 = list()
  xytab = map(lambda x : round(x,1),list(xytab))
  for x in xytab :
    activityLevels1 += [x] * 10
    activityLevels2 += xytab
  
  if shuffled :
    # keeping 0,0 as the first test and shuffling the rest
    zipper = zip(activityLevels1[1:],activityLevels2[1:])
    rnd.shuffle(zipper) # to have it with a random seed : rnd.shuffle(zipper,rnd.random)
    zipper = zip(*zipper)
    activityLevels1[1:] = zipper[0]
    activityLevels2[1:] = zipper[1]
  
  activityLevels = np.array([activityLevels1,activityLevels2]) # for both channels
  nbTimes = len(activityLevels[0])  
  CSNrate= gCSN * activityLevels + np.ones((nbTimes)) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((nbTimes)) * PTNFR[0]
  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop = {'CSN':[(),()],'PTN':[(),()]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in range(2):
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in range(2):
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))

  #-------------------------
  # log-related variables
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,nbTimes)) # for the 100 tests

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(nbRecord):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'step , '
  for i in range(nbRecord):
    for N in NUCLEI:
      frstr += N+' ('+str(i)+') , '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)
  
  #----------------------------------
  # Loop over all the steps of the test
  #----------------------------------
  for timeStep in range(nbTimes):
    print "\n\tTime step %d\n" % timeStep
    key = (round(activityLevels[0][timeStep],1),round(activityLevels[1][timeStep],1))
    #-------------------------
    # measures                                                                                                                                                  
    #-------------------------
    spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels

    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'

    for i in range(nbRecord):
      for N in NUCLEI:
        spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(timeStep)+'_'+antagStr+N, "to_file": False, 'start':offsetDuration + timeStep*(offsetDuration+simDuration),'stop':(timeStep+1)*(offsetDuration+simDuration)})
        nest.Connect(Pop[N][i], spkDetect[i][N])

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    '''
    if 'Fake' in globals():
      if 'CSN' in Fake:
        nest.SetStatus(Fake['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Fake['CSN'][1],{'rate':CSNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      if 'PTN' in Fake:
        nest.SetStatus(Fake['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Fake['PTN'][1],{'rate':PTNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    else:
      nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
      nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
      nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    '''

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    nest.Simulate(simDuration+offsetDuration)

    for i in range(nbRecord):
      print '------ Channel',i,'-------'
      for N in NUCLEI:
        #strTestPassed = 'NO!'
        expeRate[N][i,timeStep] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
        print 't('+str(timeStep)+')',N,':',expeRate[N][i,timeStep],'Hz'
        frstr += '%f , ' %(expeRate[N][i,timeStep])

    strTestPassed = 'YES!'
    choosenChannel = "0"
    # for the no input activity in both channels
    if activityLevels[0][timeStep]==activityLevels[1][timeStep]==0 :
      for i in range(params['nbCh']):
        if expeRate['GPi'][i,timeStep]<FRRNormal['GPi'][0] or expeRate['GPi'][i,timeStep]>FRRNormal['GPi'][1]:
          strTestPassed = 'NO!'
          choosenChannel = "3"
        else :
          choosenChannel = "0"
      meanRestGPi = expeRate['GPi'][:,timeStep].mean()
    # for equal & non zero input activity
    elif activityLevels[0][timeStep] == activityLevels[1][timeStep] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity higher then channel 2 input activity
    elif activityLevels[0][timeStep] > activityLevels[1][timeStep] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity lower then channel 2 input activity      
    elif activityLevels[0][timeStep] < activityLevels[1][timeStep] :
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    
    print '------ Result ------'
    print expeRate['GPi'][0,timeStep],'Hz',expeRate['GPi'][1,timeStep],'Hz',strTestPassed

    chan1Val = activityLevels[0][timeStep]
    chan2Val = activityLevels[1][timeStep]
    print "channel1 : ",chan1Val
    print "channel2 : ",chan2Val
    print "gagnant = ", choosenChannel
    print "--------------------"
    if trials_dico.has_key(key) :
      trials_dico[key].append(choosenChannel)
    else :
      trials_dico[key] = [choosenChannel]
      
    # write measured firing rates in csv file
    #frstr+='\n'
    #firingRatesFile.writelines(frstr)

    os.system("rm -f log/*.gdf")
    
    #-------------------------
    # Displays
    #-------------------------
    '''
    if showRasters and interactive:
      for i in range(nbRecord):
        displayStr = ' Channel '+str(i)
        displayStr+=' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
        #for N in NUCLEI:
        for N in ['MSN','STN']:
          #nest.raster_plot.from_device(spkDetect[i][N],hist=True,title=N+displayStr)
          nest.raster_plot.from_device(spkDetect[i][N],hist=False,title=N+displayStr)

      nest.raster_plot.show()
    '''
  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()

  return xytab,trials_dico


#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGurneyTestGenericReZero(frtrials_dico,ratio=1.5,shuffled=True,xytab=np.arange(0.,1.,0.1),showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag='',reversing=False,shutOffset=False,offset=200,simuTime=800,sameVal=None,constantChan=None,rezero=False):
  
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = float(offset)
  simDuration = float(simuTime) # ms
  
  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  loadLG14params(params['LG14modelID'])

  # VARIABLES :
  selected = []
  steps = []
  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nepa in necessaryParams:
    if nepa not in params:
      print "Missing parameter:",nepa
      exit()

  nbRecord = 2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels
      
  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()
  connectBG_MC(antagInjectionSite,antag)

  #-------------------------
  # prepare the firing rates of the inputs for the 5 steps of the experiment
  #-------------------------  
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  
  activityLevels1 = list()
  activityLevels2 = list()
  xytab = map(lambda x : round(x,1),list(xytab))
  for x in xytab :
    # to associate every nb with every nb
    activityLevels1 += [x] * len(xytab)
    activityLevels2 += xytab
  
  if shuffled :
    # keeping 0,0 as the first test and shuffling the rest
    zipper = zip(activityLevels1[1:],activityLevels2[1:])
    rnd.shuffle(zipper) # to have it with a random seed : rnd.shuffle(zipper,rnd.random)
    activityLevels = zip(*zipper)
  else :
    activityLevels = ([activityLevels1,activityLevels2]) # for both channels
  # to generate always the same tuple of activity for chan1 and 2 during length time as : sameVal=None or (chan1 activity,chan2 activity,length)
  if not sameVal is None :
    activityLevels[0] = [sameVal[0]] * sameVal[2]
    activityLevels[1] = [sameVal[1]] * sameVal[2]
  # to generate always one channel constant as constantChan=None or (chan nb,chan activity)
  if not constantChan is None :
    if sameVal :
      activityLevels[constantChan[0]] = [constantChan[1]] * sameVal[2]
    else :
      activityLevels[constantChan[0]] = [constantChan[1]] * len(xytab)
      # the other chan is a simple xytab
      activityLevels[((constantChan[0]+1)%2)] = xytab

  if activityLevels[0][0]!=0. or activityLevels[1][0]!=0. :
    activityLevels[0].insert(0,0.)
    activityLevels[1].insert(0,0.)
  nbTimes = len(activityLevels[0])
  # for generating a rest state, a (0,0) Cortical input
  if rezero :
    nbTimes *= 2
    # inserting 0.,0.
    actLevs = []
    for i in range(nbTimes) :
      if i%2==1 :
        actLevs.append((activityLevels[0][i/2],activityLevels[1][i/2]))
      else :
        actLevs.append((0.,0.))
    activityLevels = zip(*actLevs)
    
  activityLevels = np.array(activityLevels)
  CSNrate= gCSN * activityLevels + np.ones((nbTimes)) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((nbTimes)) * PTNFR[0]
  
  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  chanRange = range(2)
  if reversing :
    chanRange = reversed(chanRange)
  ActPop = {'CSN':[(),()],'PTN':[(),()]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in chanRange:
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in chanRange:
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in chanRange:
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in chanRange :
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in chanRange:
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in chanRange :
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))

 
  #-------------------------
  # log-related variables
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,nbTimes)) # for the 100 tests

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(nbRecord):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'step , '
  for i in range(nbRecord):
    for N in NUCLEI:
      frstr += N+' ('+str(i)+') , '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)
  
  #----------------------------------
  # Loop over all the steps of the test
  #----------------------------------
  for timeStep in range(nbTimes):
    steps.append(timeStep)
    print "\n\tTime step %d\n" % timeStep
    key = (round(activityLevels[0][timeStep],1),round(activityLevels[1][timeStep],1))
    
    #-------------------------
    # measures                                                                                                                                                  
    #-------------------------
    spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels

    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'

    for i in range(nbRecord):
      for N in NUCLEI:
        if shutOffset and timeStep!=0:
          start = timeStep*simDuration 
          stop = (timeStep+1)*simDuration
        else :
          start = offsetDuration + timeStep*(offsetDuration+simDuration) 
          stop = (timeStep+1)*(offsetDuration+simDuration)

        spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(timeStep)+'_'+antagStr+N, "to_file": False, 'start': start ,'stop':stop})
        nest.Connect(Pop[N][i], spkDetect[i][N])

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    if shutOffset and timeStep != 0:
      nest.Simulate(simDuration)
    else :
      nest.Simulate(simDuration+offsetDuration)

    for i in range(nbRecord):
      print '------ Channel',i,'-------'
      for N in NUCLEI:
        #strTestPassed = 'NO!'
        expeRate[N][i,timeStep] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
        print 't('+str(timeStep)+')',N,':',expeRate[N][i,timeStep],'Hz'
        frstr += '%f , ' %(expeRate[N][i,timeStep])

    strTestPassed = 'YES!'
    choosenChannel = "0"
    # for the no input activity in both channels
    if key[0]==key[1]==0 :
      for i in range(params['nbCh']):
        if expeRate['GPi'][i,timeStep]<FRRNormal['GPi'][0] or expeRate['GPi'][i,timeStep]>FRRNormal['GPi'][1]:
          strTestPassed = 'NO!'
          choosenChannel = "3"
        else :
          choosenChannel = "0"
      meanRestGPi = expeRate['GPi'][:,timeStep].mean()
    # for equal & non zero input activity
    elif key[0]==key[1] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity higher then channel 2 input activity
    elif key[0] > key[1] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity lower then channel 2 input activity      
    elif key[0] < key[1] :
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    
    print '------ Result ------'
    print expeRate['GPi'][0,timeStep],'Hz',expeRate['GPi'][1,timeStep],'Hz',strTestPassed

    chan1Val = activityLevels[0][timeStep]
    chan2Val = activityLevels[1][timeStep]
    print "channel1 : ",chan1Val
    print "channel2 : ",chan2Val
    print "gagnant = ", choosenChannel
    print "--------------------"
    selected.append(choosenChannel)
    
    if frtrials_dico.has_key(key) :
      frtrials_dico[key].append((timeStep,expeRate['GPi'][0,timeStep],expeRate['GPi'][1,timeStep],choosenChannel))
    else :
      frtrials_dico[key] = [(timeStep,expeRate['GPi'][0,timeStep],expeRate['GPi'][1,timeStep],choosenChannel)]
      
    # write measured firing rates in csv file
    #frstr+='\n'
    #firingRatesFile.writelines(frstr)

    os.system("rm -f log/*.gdf")
    
    #-------------------------
    # Displays
    #-------------------------
  # END FOR
    
  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()
  
  return steps,expeRate['GPi'],selected,activityLevels,frtrials_dico

#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
# TODO To combine with the previous one and remove
def checkGurneyTestGenericZero(frtrials_dico,ratio=1.5,shuffled=True,xytab=np.arange(0.,1.,0.1),showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag=''):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = 200.
  simDuration = 800. # ms
  loadLG14params(params['LG14modelID'])

  # VARIABLES :
  selected = []
  steps = []
  # To remember the activity at zero input activity
  zero_activity_state = dict.fromkeys(NUCLEI,{})
  
  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nepa in necessaryParams:
    if nepa not in params:
      print "Missing parameter:",nepa
      exit()

  nbRecord = 2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels
      
  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()
  connectBG_MC(antagInjectionSite,antag)
        
  #-------------------------
  # prepare the firing rates of the inputs for the 5 steps of the experiment
  #-------------------------  
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  activityLevels1 = list()
  activityLevels2 = list()
  if xytab[0] != 0. :
    xytab.insert(0,0.)
  xytab = map(lambda x : round(x,1),list(xytab))
  for x in xytab :
    # to associate every nb with every nb
    activityLevels1 += [x] * len(xytab)
    activityLevels2 += xytab
  if shuffled :
    # keeping 0,0 as the first test and shuffling the rest
    zipper = zip(activityLevels1[1:],activityLevels2[1:])
    rnd.shuffle(zipper) # to have it with a random seed : rnd.shuffle(zipper,rnd.random)
    zipper = zip(*zipper)
    activityLevels1[1:] = zipper[0]
    activityLevels2[1:] = zipper[1]  
  activityLevels = np.array([activityLevels1,activityLevels2]) # for both channels
  nbTimes = len(activityLevels[0])
  CSNrate= gCSN * activityLevels + np.ones((nbTimes)) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((nbTimes)) * PTNFR[0]
  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop = {'CSN':[(),()],'PTN':[(),()]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in range(2):
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in range(2):
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  #-------------------------
  # log-related variables
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,nbTimes)) # for the 100 tests

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(nbRecord):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'step , '
  for i in range(nbRecord):
    for N in NUCLEI:
      frstr += N+' ('+str(i)+') , '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)
  
  #----------------------------------
  # Loop over all the steps of the test
  #----------------------------------
  for timeStep in range(nbTimes):
    steps.append(timeStep)
    print "\n\tTime step %d\n" % timeStep
    key = (round(activityLevels[0][timeStep],1),round(activityLevels[1][timeStep],1))
    
    # reinitialization of neuron state to initial no activity
    if timeStep!=0 :
      for i in range(nbRecord):
        for N in NUCLEI:
          print "TUPLE de taille : ",len(zero_activity_state[N][i])
          for nID,status in zero_activity_state[N][i].items() :
            print "\n------------reactualizing the zero state : \n\tneuron id=",nID," \n\tdict status=",status,"--------------\n"
            state = {'t_spike' : status['t_spike'],
                     'V_m' : status['V_m'],
                      'Ca' : status['Ca']}
            #nest.SetStatus([nID],state)
    
    #-------------------------
    # measures                                                                                                                                                  
    #-------------------------
    spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels

    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'

    for i in range(nbRecord):
      for N in NUCLEI:
        spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(timeStep)+'_'+antagStr+N, "to_file": False, 'start':offsetDuration + timeStep*(offsetDuration+simDuration),'stop':(timeStep+1)*(offsetDuration+simDuration)})
        nest.Connect(Pop[N][i], spkDetect[i][N])

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    nest.Simulate(simDuration+offsetDuration)

    for i in range(nbRecord):
      print '------ Channel',i,'-------'
      for N in NUCLEI:
        #strTestPassed = 'NO!'
        expeRate[N][i,timeStep] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
        print 't('+str(timeStep)+')',N,':',expeRate[N][i,timeStep],'Hz'
        frstr += '%f , ' %(expeRate[N][i,timeStep])
        if timeStep==0 :
          zero_activity_state[N][i] = {}
          for neuron in Pop[N][i] :
            zero_activity_state[N][i][neuron] = nest.GetStatus([neuron])[0]
            print "\n------------reactualizing the zero state : neuron = ",neuron," DATA = ", zero_activity_state[N][i][neuron],"--------------"
        print
    exit()
    strTestPassed = 'YES!'
    choosenChannel = "0"
    # for the no input activity in both channels
    if key[0]==key[1]==0 :
      for i in range(params['nbCh']):
        if expeRate['GPi'][i,timeStep]<FRRNormal['GPi'][0] or expeRate['GPi'][i,timeStep]>FRRNormal['GPi'][1]:
          strTestPassed = 'NO!'
          choosenChannel = "3"
        else :
          choosenChannel = "0"
      meanRestGPi = expeRate['GPi'][:,timeStep].mean()
    # for equal & non zero input activity
    elif key[0]==key[1] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity higher then channel 2 input activity
    elif key[0] > key[1] :
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    # for channel 1 input activity lower then channel 2 input activity      
    elif key[0] < key[1] :
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
        choosenChannel = "3"
      # both ok checking the chosen one
      elif expeRate['GPi'][0,timeStep] > expeRate['GPi'][1,timeStep] * ratio:
        choosenChannel = "2"
      elif expeRate['GPi'][1,timeStep] > expeRate['GPi'][0,timeStep] * ratio:
        choosenChannel = "1"
      else :
        choosenChannel = "0"
    
    print '------ Result ------'
    print expeRate['GPi'][0,timeStep],'Hz',expeRate['GPi'][1,timeStep],'Hz',strTestPassed

    chan1Val = activityLevels[0][timeStep]
    chan2Val = activityLevels[1][timeStep]
    print "channel1 : ",chan1Val
    print "channel2 : ",chan2Val
    print "gagnant = ", choosenChannel
    print "--------------------"
    selected.append(choosenChannel)
    
    if frtrials_dico.has_key(key) :
      frtrials_dico[key].append((expeRate['GPi'][0,timeStep],expeRate['GPi'][1,timeStep]))
    else :
      frtrials_dico[key] = [(expeRate['GPi'][0,timeStep],expeRate['GPi'][1,timeStep])]
      
    # write measured firing rates in csv file
    #frstr+='\n'
    #firingRatesFile.writelines(frstr)

    os.system("rm -f log/*.gdf")
    
    #-------------------------
    # Displays
    #-------------------------
  # END FOR
    
  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()
  
  return steps,expeRate['GPi'],selected,activityLevels,frtrials_dico

#-----------------------------------------------------------------------
def main():
  rasters = False
  WITH_GDF = False
  
  if len(sys.argv) >= 2:
    print "Command Line Parameters"
    paramKeys = ['LG14modelID',
                 'nbMSN',
                 'nbFSI',
                 'nbSTN',
                 'nbGPe',
                 'nbGPi',
                 'nbCSN',
                 'nbPTN',
                 'nbCMPf',
                 'GMSN',
                 'GFSI',
                 'GSTN',
                 'GGPe',
                 'GGPi', 
                 'IeGPe',
                 'IeGPi',
                 'inDegCSNMSN',
                 'inDegPTNMSN',
                 'inDegCMPfMSN',
                 'inDegFSIMSN',
                 'inDegMSNMSN', 
                 'inDegCSNFSI',
                 'inDegPTNFSI',
                 'inDegSTNFSI',
                 'inDegGPeFSI',
                 'inDegCMPfFSI',
                 'inDegFSIFSI',
                 'inDegPTNSTN',
                 'inDegCMPfSTN',
                 'inDegGPeSTN',
                 'inDegCMPfGPe',
                 'inDegSTNGPe',
                 'inDegMSNGPe',
                 'inDegGPeGPe',
                 'inDegMSNGPi',
                 'inDegSTNGPi',
                 'inDegGPeGPi',
                 'inDegCMPfGPi',
                 ]
    if len(sys.argv) == len(paramKeys)+1:
      print "Using command line parameters"
      print sys.argv
      i = 0
      for k in paramKeys:
        i+=1
        params[k] = float(sys.argv[i])
    else :
      print "Incorrect number of parameters:",len(sys.argv),"-",len(paramKeys),"expected"

  nest.set_verbosity("M_WARNING")
  
  #execTime = time.localtime()
  #timeStr = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])+':'+str(execTime[5])
  os.system("rm -rf log/*")
  
  score = np.zeros((2))
  
  
  score += checkAvgFR_MC(params=params,antagInjectionSite='none',antag='',showRasters=True)
  '''
  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
    score += checkAvgFR_MC(params=params,antagInjectionSite='GPe',antag=a)

  for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
    score += checkAvgFR_MC(params=params,antagInjectionSite='GPi',antag=a)
  '''

  #score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=0.2,PActivePTN=0.2)
  #score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=1.,PActivePTN=1.)

  #-------------------------
  print "******************"
  print "* Score:",score[0],'/',score[1]
  print "******************"

  #-------------------------
  # log the results in a file
  #-------------------------
  #res = open('log/OutSummary_'+timeStr+'.txt','a')
  res = open('log/OutSummary.txt','a')
  for k,v in params.iteritems():
    res.writelines(k+' , '+str(v)+'\n')
  res.writelines("Score: "+str(score[0])+' , '+str(score[1]))
  res.close()

  res = open('score.txt','w')
  res.writelines(str(score[0])+'\n')
  res.close()
  
  # Combining in one and removing the numerous number of GDF files
  if WITH_GDF :
    os.system("sh ./gdf_concat.sh")
  else : # removing every generated gdf file
    os.system("rm -f log/*.gdf")

#---------------------------
if __name__ == '__main__':
  main()
