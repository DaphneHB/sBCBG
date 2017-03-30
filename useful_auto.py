# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:29:18 2017

@author: daphnehb
"""

from testFullBG import *
import os
import plot_tools as pltT
import io_gest as io

def generate_table(model) :
  params['LG14modelID'] = model
  score = np.zeros((2))
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
  
  print "Plot table generation for the model ",str(model)
  
  global inDegree_boarders  
  pltT.plot_inDegrees_boarders_table(inDegree_boarders,model)
  
def generate_margin_plot(path=None, model=None, antag=None, norm=False, score=0, limit=-1):
  if path is None :
    path = os.getcwd()
  filename = 'allFiringRates'
  io.concat_data(outFile=filename,dataPath=path,model=model,score=score,limit=limit)
  pltT.plot_margins_and_simus(filename=filename,antag=antag,norm=norm,model=model)

'''
Generate the plot(s) for analyzing the Ie and G params on nucleus
variables must be a list of Ie and/or G
nucleus must be a list of BG nucleus
each variables is then associated with each nucleus
a tuple (variable,nucleus) define a figure
'''    
def generate_param_analyze(variables, nucleus, path=None, score=0, model=None, separated=False) :
  if path is None :
    path = os.getcwd()
  # whether we want to plot for multiple params
  if len(variables) == len(nucleus) == 1 :
    pltT.plot_score_ratio(variables[0], nucleus[0], dataPath=path, score=score, model=model, axis=None)
  elif not separated  :
    nbN = len(nucleus)
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
        pltT.plot_score_ratio(variables[j], nucleus[i], dataPath=path, score=score, model=model, axis=ax)
  else :
    nbN = len(nucleus)
    nbV = len(variables)
    for i in range(nbN) :
      for j in range(nbV) :
        pltT.plot_score_ratio(variables[j], nucleus[i], dataPath=path, score=score, model=model, axis=None)
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
def generate_fr_by_param(N, V, interv) :
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
    score,score_max = checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=False)
    # getting the FR in firingRates.csv
    with open("log/firingRates.csv") as frFile :
      allFRdata = frFile.readline().rstrip().split(",") # no antag => one line
    # Retrieving the FR of the nucleus N
    fr = float(allFRdata[3:-1][NUCLEI.index(N)])
    vals.append((val,fr,score))
        
  pltT.plot_fr_by_var(n_var, vals, score_max)
  
#generate_table(2)
#generate_margin_plot(path="/home/daphnehb/OIST/SangoTests/model2/copyBG",limit=50,score=3)
#generate_param_analyze(['Ie', 'G'], ['GPi','GPe','FSI'], path="/home/daphnehb/OIST/SangoTests/model2/copyBG",separated=True)

generate_fr_by_param('GPi','G', (0.5,26,0.1))
'''
for N in NUCLEI :
  # for each nucleus generating G plot
  generate_fr_by_param(N,'G', (0.5,26,0.1))
  if "GP" in N :
    generate_fr_by_param(N,'Ie', (5,30,1))
'''