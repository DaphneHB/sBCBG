# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:29:18 2017

@author: daphnehb
"""

from testFullBG import *

import os
import re
import plot_tools as pltT
import io_gest as io
from datetime import datetime

'''
To download the accurate parametrization for the current model
'''
def switch_model(model) :
  if model==0 :
    from modelParams0 import params
  elif model==1 :
    from modelParams1 import params
  elif model==2 :
    from modelParams2 import params
  elif model==3 :
    from modelParams3 import params
  elif model==4 :
    from modelParams4 import params
  elif model==5 :
    from modelParams5 import params
  elif model==6 :
    from modelParams6 import params
  elif model==7 :
    from modelParams7 import params
  elif model==8 :
    from modelParams8 import params
  elif model==9 :
    from modelParams9 import params
  elif model==10 :
    from modelParams10 import params
  elif model==11 :
    from modelParams11 import params
  elif model==12 :
    from modelParams12 import params
  elif model==13 :
    from modelParams13 import params
  elif model==14 :
    from modelParams14 import params

def launching_exec_by_models(validFile,models=np.arange(0,15,1),score=0) :
  # removing the previous tests
  os.system("rm -rf "+validFile)
  # TODO charge paramFIlePath modelParams file if it exist
  for mod in models :
    params['LG14modelID'] = mod
    print "Generating for model #" + str(mod)
    os.system("rm -r log/*")
    scoreTab = np.zeros((2))
    scoreTab += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
    if with_antag :      
      for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
        scoreTab += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)
    
      for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
        scoreTab += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)
    if scoreTab[0] < score :
      continue
      
def launching_exec_by_intervalle(n_var,interv,validFile,model=0,with_antag=False) :
  os.system("rm -rf " + validFile)
  # TODO charge modelParams file for params
  params['LG14modelID'] = model
  for val in np.arange(*interv) :
    params[n_var] = float(val)
    print "****** Simulation for " + n_var + " = " + str(val)
    # emptying the log file before each simulation
    os.system("rm -r log/*")
    # launching simulations
    checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
    if with_antag :      
      for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
        checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)
    
      for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
        checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)
    
'''
Compute the inDegree for a certain model to be in the ranges the relative way #9 is

According to this formula
For model #9
  - for each connection k:
      the Liénard model defines inDegree(#9)^{max}_{k} and
inDegree(#9)^{min}_{k}
      I chose some inDegree(#9)_k value
For any other model #N, set inDegree(#N)_{k} the following way:
  inDegree(#N)_{k} = inDegree(#N)^{min}_k + ( inDegree(#9)_k -
inDegree(#9)^{min}_k ) * (inDegree(#N)^{max}_k - inDegree(#N)^{min}_k) /
(inDegree(#9)^{max}_k - inDegree(#9)^{min}_k)
'''
def compute_inDegree(model) :
  from modelParams9 import params as params9
  print params9
  exit
  if model==9 :
    return pltT.retreive_inDegree(model,prms=params9)
  else :
    nine_boarders = pltT.retreive_inDegree(9,prms=params9)
    model_boarders = pltT.retreive_inDegree(model)
    for nameTgt in NUCLEI :
      for nameSrc in nbSim.keys() :
        key = nameSrc + "->" + nameTgt
        if model_boarders.has_key(key) and nine_boarders.has_key(key) :
          # according to the formula :
          nine_min,nine_max,nine_val = nine_boarders[key]
          model_min,model_max,model_val = model_boarders[key]
          new_inDegree = model_min + (nine_val - nine_min) * (model_max - model_min) / (nine_max - nine_min)
          # saving        
          model_boarders[key] = model_min,model_max, new_inDegree
    return model_boarders
        

'''
Retrieve the parametrization used
Stocked in the global variable params or modelParams.py file
'''
def get_params(paramsFromFile=None,remove=[]) :
  if paramsFromFile is None :
    # getting params from params dict
    global NUCLEI, params
    legend = []
    for N in NUCLEI :
      # getting the gain for this Nucleus
      prm = "G" + N
      if not prm in remove :
        val = prm + " = " + str(params[prm])
        legend.append(val)
      # also getting the input current
      if N=="GPe" or N=="GPi" :
        prm = "Ie" + N
        if not prm in remove :
          val = prm + " = " + str(params[prm])
          legend.append(val)
    return legend
  else :
    return io.get_param_from_file(paramsFromFile)
    
'''
Generate the inDegree table from each nucleus to each nucleus for the given model number
show the plot
'''
def generate_table(model,save=False) :
  switch_model(model)
  params['LG14modelID'] = model
  print "Generating inDegree Table for model " + str(model)
  
  inDegree_boarders = pltT.retreive_inDegree(model)
  
  print "Plot table generation for the model ",str(model)
  
  
  if save :
    filename = "tables/inDegreeTable#" + str(model) + ".png"
  else :
    filename = None
  pltT.plot_inDegrees_boarders_table(inDegree_boarders,model=model,filename=filename)
  

'''
Generate the plot with the margin to analyze the in-range data and their scores
Useful to see if the data are in the high-range or low-range

glob defines whether or not we analyse multiple data or only one firingRates.csv file
if glob is true : the path describe the path to the firingRates.csv simple file
path is the path where the to be concatenate data can be found
model is the model to analyse, none: for every simulated models
antag defines whether or not we want to analyse the antagonist injection as well
norm (not with antag) : to translate the margin ranges into 0-1
score : the min expected obtained score to analyze
separated : plotting separatedly or in the same figure
'''
def generate_margin_plot(glob=True,path=None, model=None, antag=None, norm=False, score=0, limit=-1,separated=False):
  if path is None :
    path = os.getcwd()
  if antag=='none' :
    antag = None
  filename = 'allFiringRates'
  if glob :
    io.concat_data(outFile=filename,dataPath=path,model=model,score=score,limit=limit)
  else : # coping the input firingRates.csv file to log as allFiringRates.csv
    os.system("cp " + os.path.join(path,"log/firingRates.csv") + " " + os.getcwd() +"/log/allFiringRates.csv")
  pltT.plot_margins_and_simus(filename=filename,antag=antag,norm=norm,model=model,separated=separated)


'''
Generate the plot(s) for analyzing the Ie and G params on nucleus
variables must be a list of Ie and/or G
nucleus must be a list of BG nucleus
each variables is then associated with each nucleus
a tuple (variable,nucleus) define a figure
'''    
def generate_param_score_analyze(variables, nuclei, path=None, score=0, model=None, separated=True, save=False) :
  if path is None :
    path = os.getcwd()
  # whether we want to plot for multiple params
  if len(variables) == len(nuclei) == 1 :
    plotName = "plots/scoreRatio" + variables[0] + nuclei[0] + "#" + str(model) + ".png"
    pltT.plot_score_ratio(variables[0], nuclei[0], dataPath=path, score=score, model=model, axis=None,save=plotName)
  elif not separated  :
    nbN = len(nuclei)
    nbV = len(variables)
    fig, axes = pltT.plt.subplots(nrows=nbN , ncols=nbV)
    fig.canvas.set_window_title("Scores of simulations under param values")
    for i in range(nbN) :
      for j in range(nbV) :
        ax = axes[i]
        if type(ax) is np.ndarray :
          ax = ax[j]
        if i==j==0 :
          ax.set_ylabel('Number of simulations')
        pltT.plot_score_ratio(variables[j], nuclei[i], dataPath=path, score=score, model=model, axis=ax,save=None)
    if save :
      nucleiStr = "_" + "+".join(nuclei) + "_" + "+".join(variables)
      fig.savefig("plots/scoreRatio" + nucleiStr + "#" + str(model) + ".png")
  else :
    nbN = len(nuclei)
    nbV = len(variables)
    for i in range(nbN) :
      for j in range(nbV) :
        plotName = "plots/scoreRatio" + variables[j] + nuclei[i] + "#" + str(model) + ".png"
        pltT.plot_score_ratio(variables[j], nuclei[i], dataPath=path, score=score, model=model, axis=None,save=plotName)
        
  if not save :
    pltT.plt.show()
    

'''
Generate for a speccific variable v Ie or G 
for a specific nucleus n the FR fct of the variable
when the variable vary in (min,max) intervalle
The variation is with a certain step
N is the string corresponding to a BG nucleus
V is the variable Ie or G we want to make vary (string)
interv is the expected details of the wanted intervalle (tuple : (min,max,step))
'''
def generate_fr_by_param(N, V, interv, model=0, with_antag=False) :
  global NUCLEI
  if not N in NUCLEI :
    print "------------ ERROR : Nucleus " + N + " does not exist"
    return 1
  if not V=='G' and not V=='Ie' :
    print "------------ ERROR : Variable " + V + " does not exist"
    return 1    
  if len(interv) != 3 :
    print "------------ ERROR : Wrong Intervalle"
    return 1
  if not (type(N) is str and type(V) is str) :
    print "------------ ERROR : Nucleus and Parameter must be string args"
    return 1
  n_var = V + N
  if not params.has_key(n_var) :
    print "------------ ERROR : Parameter " + n_var + " does not exist"
    return 1
    
  params['LG14modelID'] = model
  print "Simulating under " + V + " for " + N
  vals = []
  for val in np.arange(*interv) :
    params[n_var] = float(val)
    print "****** Simulation for " + n_var + " = " + str(val)
    # emptying the log file before each simulation
    os.system("rm -r log/*")
    score = np.zeros((2))
    # launching simulations
    
    score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
    if with_antag :      
      for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
        score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)
    
      for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
        score += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)
        
    sc, score_max = score
    # getting the FR in firingRates.csv
    with open("log/firingRates.csv") as frFile :
      allFRdata = frFile.readline().rstrip().split(",") # no antag => one line
    # Retrieving the FR of the nucleus N
    fr = float(allFRdata[3:-1][NUCLEI.index(N)])
    vals.append((val,fr,sc))
  
  pltT.plot_fr_by_var(n_var, vals, score_max,interv, model)
  

'''
Analyze mutually param1, param2, param3
'''
def generate_param_analyze(param1,param2,param3=None,save=False,path=os.getcwd(), score=0,model=None) :
  pltT.plot_param_by_param(param1,param2,param3=param3,save=save,dataPath=path, score=score,model=model)
  
  
  
def generate_models_ranges_tab(parametrization=None,to_generate=True, validationPath=os.getcwd(),paramFilePath=None, models=np.arange(0,15,1), score=0, with_antag=False,save=False) :
  allFRdata = {}
  validFile = os.path.join(validationPath,"validationArray.csv")
    
  if to_generate :
    launching_exec_by_models(validFile,models=models,score=score,with_antag=with_antag)
  #path = os.path.dirname(paramFilePath)
  with open(validFile, 'r') as frFile :
    FRdata = frFile.readlines()  # list of lines
  nmodels = []
  for mod in models :
    try:
      data = filter(lambda x : ("#"+str(mod)+" ") in x, FRdata)
    except  IndexError :
      continue
    nmodels.append(mod)
    allFRdata[mod] = data
  legend = get_params(paramsFromFile=paramFilePath)  
  print "Plot array simu generation for models ",str(nmodels)
  # plotting with plot_tools file
  filename = None
  if save :
      filename = "plots/" + str(datetime.now()) + "modelsValidation.png"
  pltT.plot_models_ranges(allFRdata, legend, filename=filename, models=nmodels)
  
  
 
def generate_gap_from_range_local(n_var, interv,save=False,to_generate=True,model=0,paramFilePath=None,pathToFile=os.getcwd(),removing=[],with_antag=False) :
  global NUCLEI
  validFile = os.path.join(pathToFile,"validationArray.csv")
  print "Simulating for " + n_var
  if to_generate :
    if len(interv) != 3 :
      print "------------ ERROR : Wrong Intervalle"
      return 1
    print "****************Simulation with Model %d******************" % model
    if not params.has_key(n_var) :
      print "------------ ERROR : Parameter " + n_var + " does not exist"
      return 1
    launching_exec_by_intervalle(n_var,interv,validFile,model=model,with_antag=with_antag)
  # getting the FR gaps in validationArray file
  vals = np.arange(*interv)
  results = io.read_validationArray_values(pathToFile=pathToFile,model=model,with_antag=with_antag)
  removing.append(n_var)
  print "Plot array simu generation for the model ",str(model)
  legend = get_params(paramsFromFile=paramFilePath,remove=removing)
  # plotting with plot_tools file
  filename = None
  if save :
      filename = "plots/" + str(datetime.now()) + "gapPlot_"+n_var+"_model"+str(model)+".png"
  pltT.plot_gap_from_range(vals, n_var, interv, results, model, filename=filename,param=legend)


# completely TODO
def generate_gap_from_range_global(n_var,save=False,model=0,paramFilePath=None,pathToFile=os.getcwd(),removing=[]) :
  global NUCLEI
  validFile = os.path.join(pathToFile,"validationArray.csv")
  print "Simulating for " + n_var
  print "****************Simulation with Model %d******************" % model
  if not params.has_key(n_var) :
    print "------------ ERROR : Parameter " + n_var + " does not exist"
    return 1
  pltT.get_data
  launching_exec_by_intervalle(n_var,interv,validFile,model=model)
  # getting the FR gaps in validationArray file
  vals = np.arange(*interv)
  results = io.read_validationArray_values(pathToFile=pathToFile,model=model)
  removing.append(n_var)
  print "Plot array simu generation for the model ",str(model)
  legend = get_params(paramsFromFile=paramFilePath,remove=removing)
  # plotting with plot_tools file
  filename = None
  if save :
      filename = "plots/" + str(datetime.now()) + "gapPlot_"+n_var+"_model"+str(model)+".png"
  pltT.plot_gap_from_range(vals, n_var, interv, results, model, filename=filename,param=legend)
  
  
  
def generate_best_score_comp(figAxes=None,path=os.getcwd(), score=0, model=None, separated=True, save=False) :
  global NUCLEI
  # whether we want to plot for multiple params
  nbN = len(NUCLEI)
  parameters = ["G" + N for N in NUCLEI] + ["IeGPe","IeGPi"]
  # retrieving the nucleus score per simulation according to the value and param
  paramData = pltT.get_data_by_model(parameters,model=model,path=path,score=score)
  if paramData == {}:
    print "----------------- ERREUR : No data to analyse"
    return 1
  simu_color = {}  # simu_color = {simu : color}
  if not separated  :
    if figAxes is None :
      fig, axes = pltT.plt.subplots(nrows=1 , ncols=len(parameters))
      fig.canvas.set_window_title("Best scores for every simulations #" + str(model))
    else :
      fig,axes = figAxes
      axes.set_title("#" + str(model))
    for i,p in enumerate(parameters) :
      ax = axes[i]
      if i==0 :
        ax.set_ylabel('Score')
      simu_color = pltT.plot_score_by_value(p, paramData[p], simu_color, model=model, axis=ax,filename=None)
    if save :
      nucleiStr = "_" + "+".join(parameters)
      fig.savefig("plots/" + str(datetime.now()) + "bestScores" + nucleiStr + "#" + str(model) + ".png")
  else :
    for i,p in enumerate(parameters) :
      plotName = "plots/" + str(datetime.now()) + "bestScores" + p + "#" + str(model) + ".png"
      simu_color = pltT.plot_score_by_value(p, paramData[p], simu_color, model=model, axis=None,filename=plotName)
  if not save :
    pltT.plt.show()

'''
Generate the multi pie chart plot for Gurney&Prescott test
Test if the right channel is selected for the action according to the input
Plus, write a file with the results (filename = "dualchanCompetition.csv")

Arguments :
NbTrials is the number of trials (1 trial = 100 simulations)
pathToData is to know where to save the file
model is model number to test
antag is none/all/GPe_AMPA/GPe_GABAA/GPe_NMDA/GPe_AMPA+GABAA/GPi_NMDA/GPi_AMPA/GPi_GABAA/GPI_NMDA+AMPA
'''
def generate_GurneyTest(generate=True,filename=None,NbTrials=5,pathToData=os.getcwd(),model=0,antag='none') :
  from testChannelBG import *
  
  trials_dico = {}
  xytab = None
  if generate or filename is None:
    for trial in range(NbTrials) :
      os.system("rm -r log/*")
      print "############################## TRIAL no. %d ###########################" % trial
      # launch GurneyTestsgeneric
      xytab,this_trial = checkGurneyTestGeneric(trials_dico,showRasters=True,params=params,PActiveCSN=0.2,PActivePTN=0.2)
      trials_dico.update(this_trial)
    if filename is None :
      filename = str(datetime.now()) + "dualchanCompetition.csv"
    io.write__2chan_file(xytab,trials_dico,filename,pathToFile=pathToData,model=model,antag=antag)
  else :
    xytab,trials_dico = io.read_2chan_file(filename,pathToFile=pathToData,model=model,antag=antag)
  print trials_dico
  trials_dico = pltT.dualchanFileToPercentages(trials_dico)
  print trials_dico
  pltT.plot_multichan_pieChart(xytab,trials_dico)

'''
for md in range(0,15) :
  generate_table(md,save=True)
'''

#generate_table(2)
#generate_margin_plot(glob=False,antag='none',path="/home/daphnehb/OIST/sBCBG3/",limit=100,score=0)
#generate_margin_plot(glob=True,antag='all',path="/home/daphnehb/OIST/SangoTests/model5/2017_4_13/",limit=5,score=0)
#generate_param_score_analyze(['G','Ie'], ['MSN','FSI','GPe','GPi','STN'],score=0, save=False,separated=True,path="/home/daphnehb/OIST/SangoTests/model1/2017_4_21")
#generate_param_analyze("GMSN","GFSI",param3="GSTN", score=11,save=False,path="/home/daphnehb/OIST/SangoTests/model1/2017_4_21",model=1)
#generate_fr_by_param('GPe','Ie', (5.,13.,1),with_antag=True)
#generate_fr_by_param('GPi','Ie', (5.,15.,1),with_antag=True)
#generate_fr_by_param('MSN','G', (4.,7.,0.25),with_antag=True)
'''
generate_fr_by_param('GPi','G', (5.,15.,0.5),with_antag=True)
generate_fr_by_param('STN','G', (1.2,1.5,0.01),with_antag=True)
generate_fr_by_param('GPe','G', (0.01,0.5,0.01),with_antag=True)
#generate_fr_by_param('FSI','G', (0.5,1.8,0.01),with_antag=True)
'''

#generate_fr_by_param('MSN','G', (4.,6.,0.1), model=1,with_antag=True)

#params['LG14modelID'] = 3
#params['GMSN'] = 5.7
#generate_fr_by_param('GPi','G', (0.5,6.,0.1),with_antag=True)
'''
for N in NUCLEI :
  # for each nucleus generating G plot
  #generate_fr_by_param(N,'G', (0.5,26,0.1))
  if "GPi" in N :
    generate_fr_by_param(N,'Ie', (5,30,1))
'''
#params['LG14modelID'] = 3
'''
for g in np.arange(4.,6.,0.1) :
  params['GMSN'] = g
  print "GGGGGGGGGG = ",g
  generate_models_ranges_tab(parametrization=params, to_generate=True,with_antag=True,save=True)
'''
#params['GMSN'] = 4.8
#params['GFSI'] = 1.2
#params['GSTN'] = 1.33
#params['GGPe'] = 1.0
#generate_gap_from_range_local("IeGPe",(6.,7.,1.),to_generate=True,model=14,save=False,removing=[],with_antag=True)
#generate_models_ranges_tab(parametrization=params, to_generate=False,with_antag=True,save=False)
#generate_best_score_comp(figAxes=None,path="/home/daphnehb/OIST/SangoTests/model5/2017_4_13/", score=0, model=5, separated=True, save=False)

#for items in compute_inDegree(14).items() :
#  print items

generate_GurneyTest(NbTrials=1)
#generate_GurneyTest(generate=False, filename="dualchanCompetition.csv",pathToData="/home/daphnehb/OIST/sBCBG3/data/tests",model=1,antag="none")
#per = pltT.dualchanFileToPercentages(di)
#print per
#pltT.plot_multichan_pieChart(np.arange(0,1.1,0.1),per)