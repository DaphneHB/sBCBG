# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 17:44:34 2017

@author: daphnehb
"""

import re
import os
import random
import commands
import numpy as np
from operator import add

from pylab import cm

from LGneurons import NUCLEI, nbSim, FRRNormal, dataPath, FRRAnt, recType
import modelParams as mparams

'''
Display on a plot the reason why a slope or bars cant ba displayed
'''
def plot_print_wrong(ax,reason) :
  if not ax is None :
    x0, x1, y0, y1 = ax.axis()
    ax.text(x1 - x0 + 2., y1 - y0, reason, ha="center", family='sans-serif', size=12)
  return 1
  

'''
According to the norm param, the function will normalize 
from certain ranges in FRRNormal to 0-1 the input in frlist.
'''
def normalize(frlist,norm=False) :
  if (not norm or len(FRRNormal)!=len(frlist)) :
    return frlist
  new_list = [0] * 5
  global NUCLEI
  for ind,N in enumerate(NUCLEI) :
    new_list[ind] = (frlist[ind]-FRRNormal[N][0])/float(FRRNormal[N][1] - FRRNormal[N][0])
  return new_list
  

'''
Function which organise the inDegrees in a 5*nb BG Nulcei*(nb BG nuclei + other Nuclei)
table_dict : {Src->Tgt: (min inDegree,max inDegree,chosen value)}
Return a matrix [[Src name,Tgt name, min inDegree, max inDegree, chosen value]]
'''
def get_inDegree_to_plot (table_dict) :
  global NUCLEI,nbSim
  nbBGNuclei = len(NUCLEI)
  nbNuclei = len(nbSim)
  # Matrix containing every data according to the nb of connections
  clust_data = np.empty((nbBGNuclei*nbNuclei,5), dtype='object')
  # for each NUCLEI we test and write if a connection exists with one of the dict key (including PTN , CSN and CMPf) which is in nbSim dict
  for i,nameTgt in enumerate(NUCLEI) :
    for j,nameSrc in enumerate(nbSim.keys()) :
      
      clust_data[i*nbNuclei+j][0] = nameSrc
      clust_data[i*nbNuclei+j][1] = nameTgt
      key = nameSrc + "->" + nameTgt

      if (table_dict.has_key(key)):
        clust_data[i*nbNuclei+j][2] = str(table_dict[key][0] if table_dict[key][0] != -1 else "unknown")
        clust_data[i*nbNuclei+j][3] = str(table_dict[key][1])
        clust_data[i*nbNuclei+j][4] = str(table_dict[key][2])
      else :
        clust_data[i*nbNuclei+j][2] = "---"
        clust_data[i*nbNuclei+j][3] = "---"
        clust_data[i*nbNuclei+j][4] = "---"
  return clust_data
  
def get_data_from_file(filterFct, filename=None,model=None) :
  # retrieving data in the input file allFiringRates.csv
  if filename is None :
    filename = 'allFiringRates'
  allFRfile = open('log/' + filename + ".csv",'r')
  allFRdata = allFRfile.readlines()
  allFRfile.close()
  # retrieving only the simu with the choosen model if there is
  if (not model is None) :
    allFRdata = filter(filterFct ,allFRdata)
  return allFRdata
  
def get_count_for_score(n_var, dataPath=os.getcwd(), score=0, model=None,axis=None) :
  val_tab = []
  varN_values = {}   # dict {score : {val : nb}}
  model_pattern = re.compile("LG14modelID.*:\ *(\d+).")
  paramVal_pattern = re.compile(n_var + ".*:\ *(\d+[\.\d*]*).*")
  for fName in os.listdir(dataPath) :
    dirPath = os.path.join(dataPath,fName)
    if os.path.isdir(dirPath) and fName.startswith("2017") :
      try:
        with open(os.path.join(dirPath, "score.txt"),"r") as scoreFile :
          obt_score = float(scoreFile.readline().rstrip())
      except Exception :
        continue
      if obt_score < float(score) :
        continue
      # If the score is ok
      # lets check the model nb by getting the modelParams
      with open(dirPath+"/modelParams.py", 'r') as paramsFile :
        Paramsdata = paramsFile.readlines()
      # only getting the results of the expected model
      if (not model is None) :
        mod = int(model_pattern.findall(filter(lambda x : model_pattern.search(x), Paramsdata)[0])[0])
        if (mod != model) :
          continue
      # get value
      try :
        val = float(paramVal_pattern.findall(filter(lambda x : paramVal_pattern.search(x), Paramsdata)[0])[0])
        if not val in val_tab :
          val_tab.append(val)
      except IndexError: # if there were no result : the variable name is wrong
        reason = "------------- ERROR : Wrong variable name [" + n_var + "]"
        print reason
        return plot_print_wrong(axis,reason)
      # extending the nb
      if (not varN_values.has_key(obt_score)) :
        varN_values[obt_score] = {val : 1} #dict(zip([float(x) for x in range(score,15)],[0.] * (15-score)))   # initializing every possible score for this value
      elif (varN_values[obt_score].has_key(val)) :
        varN_values[obt_score][val] += 1
      else :
        varN_values[obt_score][val] = 1
  # for every score # for values of the param that are not in score, put the number to 0
  for scKeys,valDict in varN_values.items() :
    for i,val in enumerate(val_tab) :
      if not valDict.has_key(val) :
        varN_values[scKeys][val] = 0
  return val_tab, varN_values
  
def get_param_param_scores(param1, param2, param3=None, dataPath=os.getcwd(), score=0, model=None) :
  param1_vals = []
  param2_vals = []
  param3_vals = []
  score_vals = []
  eachPoint = {}  # dict of list as point coordinate : score list
  
  model_pattern = re.compile("LG14modelID.*:\ *(\d+).")
  param1Val_pattern = re.compile(str(param1) + ".*:\ *(\d+[\.\d*]*).*")
  param2Val_pattern = re.compile(str(param2) + ".*:\ *(\d+[\.\d*]*).*")
  param3Val_pattern = re.compile(str(param3) + ".*:\ *(\d+[\.\d*]*).*")
  for fName in os.listdir(dataPath) :
    dirPath = os.path.join(dataPath,fName)
    if os.path.isdir(dirPath) and fName.startswith("2017") :
      try:
        with open(os.path.join(dirPath, "score.txt"),"r") as scoreFile :
          obt_score = float(scoreFile.readline().rstrip())
      except Exception :
        continue
      if obt_score < float(score) :
        continue
      # If the score is ok
      # lets check the model nb by getting the modelParams
      with open(dirPath+"/modelParams.py", 'r') as paramsFile :
        Paramsdata = paramsFile.readlines()
        # only getting the results of the expected model
      if (not model is None) :
        mod = int(model_pattern.findall(filter(lambda x : model_pattern.search(x), Paramsdata)[0])[0])
        if (mod != model) :
          continue
      # get values
      point = []
      try :
        val1 = float(param1Val_pattern.findall(filter(lambda x : param1Val_pattern.search(x), Paramsdata)[0])[0])
        param1_vals.append(val1)
        point.append(val1)
      except IndexError: # if there were no result : the variable name is wrong
        reason = "------------- ERROR : Wrong variable name [" + str(param1) + "]"
        print reason
        return plot_print_wrong(axis,reason)
      try :
        val2 = float(param2Val_pattern.findall(filter(lambda x : param2Val_pattern.search(x), Paramsdata)[0])[0])
        param2_vals.append(val2)
        point.append(val2)
      except IndexError: # if there were no result : the variable name is wrong
        reason = "------------- ERROR : Wrong variable name [" + str(param2) + "]"
        print reason
        return plot_print_wrong(axis,reason)
      if not param3 is None :
        try :
          val3 = float(param3Val_pattern.findall(filter(lambda x : param3Val_pattern.search(x), Paramsdata)[0])[0])
          param3_vals.append(val3)
          point.append(val3)
        except IndexError: # if there were no result : the variable name is wrong
          reason = "------------- ERROR : Wrong variable name [" + str(param3) + "]"
          print reason
          return plot_print_wrong(axis,reason)
      # saving every score for each point
      point = tuple(point)
      if eachPoint.has_key(point) :
        eachPoint[point].append(obt_score)
      else :
        eachPoint[point] = [obt_score]
      score_vals.append(obt_score)
  # end for
  score_vals = np.array(score_vals)
  colmap = cm.ScalarMappable(cmap=cm.hsv)
  colmap.set_array(score_vals)
  
  return param1_vals, param2_vals, param3_vals, eachPoint, score_vals, colmap