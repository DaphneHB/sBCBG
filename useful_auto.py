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

def get_params(paramsFromFile=None) :
  if paramsFromFile is None :
    # getting params from params dict
    global NUCLEI
    legend = []
    for N in NUCLEI :
      # getting the gain for this Nucleus
      prm = "G" + N
      val = prm + " = " + str(params[prm])
      legend.append(val)
      # also getting the input current
      if N=="GPe" or N=="GPi" :
        param = "Ie" + N
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
  params['LG14modelID'] = model
  print "Generating inDegree Table for model " + str(model)
  os.system("rm -r log/*")
  score = np.zeros((2))
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
  
  print "Plot table generation for the model ",str(model)
  
  global inDegree_boarders
  if save :
    filename = "tables/inDegreeTable#" + str(model) + ".png"
  else :
    filenamae = None
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
def generate_fr_by_param(N, V, interv, with_antag=False) :
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
  
  pltT.plot_fr_by_var(n_var, vals, score_max,interv)
  

'''
Analyze mutually param1, param2, param3
'''
def generate_param_analyze(param1,param2,param3=None,save=False,path=os.getcwd(), score=0,model=None) :
  pltT.plot_param_by_param(param1,param2,param3=param3,save=save,dataPath=path, score=score,model=model)
  
  
  
def generate_models_ranges_tab(parametrization=None,to_generate=True, validationPath=os.getcwd(),paramFilePath=None, models=np.arange(0,15,1), score=0, with_antag=False,save=False) :
  allFRdata = {}
  validFile = os.path.join(validationPath,"validationArray.csv")
    
  if to_generate :
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
  
  
 
def generate_gap_from_range(n_var, interv,save=False,to_generate=True,model=0,paramFilePath=None,pathToFile=os.getcwd()) :
  global NUCLEI
  print "Simulating for " + n_var
  if to_generate :
    if len(interv) != 3 :
      print "------------ ERROR : Wrong Intervalle"
      return 1
    print "****************Simulation with Model %d******************" % model
    if not params.has_key(n_var) :
      print "------------ ERROR : Parameter " + n_var + " does not exist"
      return 1
    
    os.system("rm -rf validationArray.csv")
    # TODO charge modelParams file for params
    params['LG14modelID'] = model
    for val in np.arange(*interv) :
      params[n_var] = float(val)
      print "****** Simulation for " + n_var + " = " + str(val)
      # emptying the log file before each simulation
      os.system("rm -r log/*")
      # launching simulations
      checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
    # getting the FR in firingRates.csv
  vals = np.arange(*interv)
  # pattern to get the gap
  gapPattern = re.compile("([+|-]?\d+[\.\d*]*)")
  results = dict(zip(NUCLEI,[list() for i in range(len(NUCLEI))])) # nucleus : [gap list]
    
  with open(os.path.join(pathToFile,"validationArray.csv")) as varFile :
    allVardata = varFile.readlines()
  allVardata = filter(lambda x : ("#" + str(model)) in x, allVardata)
  for line in allVardata :# Retrieving the FR of each nucleus N
    line = line.strip().split(",")
    for ind,nucl in enumerate(NUCLEI) :
      variations = line[1:-1][ind].strip().split("=")
      if variations[1] == "OK" :
        results[nucl].append(0)
      else :
        gap = float(gapPattern.findall(variations[1])[0])
        results[nucl].append(gap)
  print "Plot array simu generation for the model ",str(model)
  legend = get_params(paramsFromFile=paramFilePath)
  # plotting with plot_tools file
  filename = None
  if save :
      filename = "plots/" + "gapPlot_"+n_var+"_model"+str(model)+".png"
  pltT.plot_gap_from_range(vals, n_var, interv, results, model, filename=filename)
  
'''
for md in range(0,15) :
  generate_table(md,save=True)
'''

#generate_table(2)
#generate_margin_plot(glob=False,antag='GPe_AMPA',path="/home/daphnehb/OIST/sBCBG/",limit=100,score=3)
#generate_margin_plot(glob=True,antag='none',path="/home/daphnehb/OIST/SangoTests/model2/2017_4_6/",limit=-1,score=0)
#generate_param_score_analyze(['G','Ie'], ['MSN','FSI','GPe','GPi','STN'],score=0, save=False,separated=True,path="/home/daphnehb/OIST/SangoTests/model2/2017_4_6")
#generate_param_analyze("GMSN","IeGPi",param3=None, score=0,save=False,path="/home/daphnehb/OIST/SangoTests/model2/2017_3_29",model=2)
#generate_fr_by_param('GPe','Ie', (5.,13.,1),with_antag=True)
#generate_fr_by_param('GPi','Ie', (5.,15.,1),with_antag=True)
#generate_fr_by_param('MSN','G', (4.,7.,0.25),with_antag=True)
'''
generate_fr_by_param('GPi','G', (5.,15.,0.5),with_antag=True)
generate_fr_by_param('STN','G', (1.2,1.5,0.01),with_antag=True)
generate_fr_by_param('GPe','G', (0.01,0.5,0.01),with_antag=True)
#generate_fr_by_param('FSI','G', (0.5,1.8,0.01),with_antag=True)
'''
#params['LG14modelID'] = 2
#generate_fr_by_param('MSN','G', (5.4,6.,0.01),with_antag=True)

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
params['LG14modelID'] = 3
for g in np.arange(4.,5.75,0.25) :
  params['GMSN'] = g
  print "GGGGGGGGGG = ",g
  generate_models_ranges_tab(parametrization=params, models=[0], to_generate=False,with_antag=False,save=False)
  break
#generate_gap_from_range("GMSN",(4.,7.,0.1),to_generate=True,model=5,save=True)