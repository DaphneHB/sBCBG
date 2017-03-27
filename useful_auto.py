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
    
#generate_table(2)
generate_margin_plot(path="/home/daphnehb/OIST/SangoTests/model2/copyBG",limit=50,score=3)