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
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)
  
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
def generate_param_analyze(variables, nucleus, path=None, score=0, model=None) :
  if path is None :
    path = os.getcwd()
  # whether we want to plot for multiple params
  if len(variables) == len(nucleus) == 1 :
    pltT.plot_score_ratio(variables[0], nucleus[0], dataPath=path, score=score, model=model, axis=None)
  else :
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
  pltT.plt.show()
    
  
  
  
#generate_table(2)
#generate_margin_plot(path="/home/daphnehb/OIST/SangoTests/model2/copyBG",limit=50,score=3)
#generate_param_analyze(['Ie', 'G'], ['GPi','GPe','FSI'], path="/home/daphnehb/OIST/SangoTests/model2/copyBG")