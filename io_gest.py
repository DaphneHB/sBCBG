# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 18:10:27 2017

@author: daphnehb
"""

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
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

write_inDegree_table("0",table,"","test_table")
'''