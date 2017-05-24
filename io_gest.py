# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 18:10:27 2017

@author: daphnehb
"""

import re
import sys
import os
from LGneurons import NUCLEI, nbSim #, inDegree_boarders

IDString = "testIO"

'''
Model is the string number of the model
table_dict is the inDegree_boarders dict we want to write
dirname is the IDString of the directory created at the beginning of the simulation
filename is optional : if not there, it will be the <dirname>_<model>.txt
'''
def write_inDegree_table(model, table_dict, dirname, filename=None) :
  # if the file isnt set yet
  if (filename is None ) :
    filename = dirname + str(model)
  filename = filename + ".txt"
  
  file_text = "Model number " + str(model) + "\n"
  file_text += "Connnection (Src->Trgt) | minInDegree | maxInDegree | choosenValue \n"
  
  print "Writing the inDegree boarders table in a file (",filename,")"
  global NUCLEI,nbSim
  
  # for each NUCLEI we test and write if a connection exists with one of the dict key (including PTN , CSN and CMPf) which is in nbSim dict
  for nameScr in nbSim.keys() :
    for nameTgt in NUCLEI :
      key = nameScr + "->" + nameTgt
      # if a connection exists
      if (table_dict.has_key(key)) :
        # lets print it
        file_text += key + " : " + str(table_dict[key]) + "\n"
      # if there is no existing connection nothing is done
  
  # opening and writing in the file
  table_file = open(filename,"w")
  table_file.writelines(file_text)
  table_file.close()

'''
ATTENTION : must be cold from the same directory as log dir parent (sibling with log dir)
Create a file outFile.csv in the current log directory from the files in dataPath
which validate the same or not parametrization according to one_param
and for one or multiple model according to model param
the exp_score param is the minimum score accepted

The output file contains firing rates of multiples simulations

By default outFile is allFiringRates.csv
and the dataPath is the current directory
'''
def concat_data(dataPath=os.getcwd(), outFile="allFiringRates", model=None, score=0, limit=-1) :
  #dataPath = os.path.join(os.getcwd(),dataPath)
  modelTitle = "every models" if model is None else "model "+str(model)
  print "+++++++++++ Concataning Firing Rates for simulations of "+modelTitle+" where score >= "+str(score)
  allFRdata = []
  # to test limit equality to continue or not
  to_continue = 0 if not limit==-1 else limit + 1
  # foreach directory in dataPath
  for fName in os.listdir(dataPath) :
    dirPath = os.path.join(dataPath,fName)
    if os.path.isdir(dirPath) and fName.startswith("2017") :
      with open("score.txt","r") as scoreFile :
        obt_score = float(scoreFile.readline().rstrip())
      if obt_score < float(score) :
        continue
      # If the score is ok
      # lets check the model nb by getting the firing rates
      with open(dirPath+"/log/firingRates.csv", 'r') as frFile :
        FRdata = frFile.readlines()
      if (not model is None) :
        FRdata = filter(lambda x : ("#" + str(model)) in x ,FRdata)
      if FRdata == [] :
        continue
      if (to_continue == limit) :
        break
      to_continue += 1
      # otherwise
      allFRdata += FRdata
  with open("log/"+outFile+".csv","w") as out :
    out.writelines(allFRdata)
  print "End of concatenation. log/"+outFile+".csv created"
        
        
def get_param_from_file(paramFilePath, with_inDegree=False) :
  global NUCLEI
  legend = []  
  with open(paramFilePath,'r') as paramFile :
    paramsData = paramFile.readlines()
  for N in NUCLEI :
    # getting the gain for this Nucleus
    param = "G" + N
    paramVal_pattern = re.compile("(." + str(param) + ".*:\ *\d+[\.\d*]*),")
    val = filter(lambda x : paramVal_pattern.search(x), paramsData)[0].replace(" ","")
    val = val.replace(",\n","")
    legend.append(val)
    # also getting the input current
    if N=="GPe" or N=="GPi" :
      param = "Ie" + N
      paramVal_pattern = re.compile("(." + str(param) + ".*:\ *\d+[\.\d*]*),.*")
      val = filter(lambda x : paramVal_pattern.search(x), paramsData)[0].replace(" ","")
      val = val.replace(",\n","")
      legend.append(val)
  # TODO also get inDegree
  return legend

'''
Retrieve a dictionary where every key will correspond to a slope to which we associate its gaps
'''
def read_validationArray_values (pathToFile=os.getcwd(),model=0,with_antag=False) :
  global NUCLEI
  if with_antag :
    limit = 14
  else :
    limit = 5
  results = dict(zip(NUCLEI,[list() for i in range(len(NUCLEI))])) # nucleus : [gap list]
  with open(os.path.join(pathToFile,"validationArray.csv")) as varFile :
    allVardata = varFile.readlines()
  allVardata = filter(lambda x : ("#" + str(model)) in x, allVardata)
  if allVardata==[] :
    print "--------- ERROR : No simulation matching this model"
    exit()
  for line in allVardata :# Retrieving the FR of each nucleus N
    line = line.strip().split(",")[1:-1]
    if len(line) != limit :
      continue
    for elt in line :
#    for ind,nucl in enumerate(NUCLEI) :
      variations = elt.strip().split("=")
      nucl = variations[0]
      if not results.has_key(nucl) :
        results[nucl] = list()
      # in any case, adding the gap result
      if variations[1] == "OK" :
        results[nucl].append(0)
      else :
        gap = float(variations[1])
        results[nucl].append(gap)
  return results
  
'''
For mutli-channel selection competition
Read the data about the channels that been chosen for this trial 
for a specific (x.y) value

File format :
#model_nb
antag_type
X ; Y ; 0, 1 or 2 (choosen channel); 0, 1 or 2; 0, 1 or 2 .... [as much as the nb of trials]
.
.
.
[100 values - 10x * 10y]
'''
def read_2chan_file(pathToFile=os.getcwd(), model=0, antag="none") :
  results = {}
  with open(os.path.join(pathToFile,"dualchanCompetition.csv")) as varFile :
    allVardata = varFile.readlines()
  if not ("#" + str(model)) in allVardata[0] :
    print "--------- ERROR : No simulation matching this model (#"+str(model)+")"
    exit()
  if not antag in allVardata[1] :
    print "--------- ERROR : No simulation matching this antagonist injection ("+antag+")"
    exit()
  chan_output_dict = {}
  # reading and sadding the values to the dictionnary
  for line in allVardata[2:] :
    line = line.split(";")
    # retrieving the coordinates (tested values)
    x,y = map(float,line[0:2])
    # getting the choosen channels + removing the last \n
    chan_output_dict[(x,y)] = map(lambda v: v.strip(),line[2:-1])
  return chan_output_dict
  
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

write_inDegree_table("0",table,"","test_table")
'''
#concat_data()
