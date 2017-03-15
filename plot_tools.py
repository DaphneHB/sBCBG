# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:32:40 2017

@author: daphnehb
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from LGneurons import NUCLEI, nbSim, interactive, FRRNormal


### FUNCTIONS

'''
To plot the table of the inDegree intervalle for a specific model
Must be called at the end of the simulation
'''
def plot_inDegrees_boarders_table(table_dict, model, filename=None) :
  global NUCLEI,nbSim
  nbBGNuclei = len(NUCLEI)
  nbNuclei = len(nbSim)
  # Matrix containing every data according to the nb of connections
  clust_data = np.empty((nbBGNuclei*nbNuclei,5), dtype='object')
  # Labels for the column of the table
  collabel = ("Src", "Target", "Min inDegree", "Max inDegree", "Choosen Value")
  
  nrows, ncols = len(clust_data)+1, len(collabel)
  hcell, wcell = 0.2, 1.
  hpad, wpad = 0, 1.5

  fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
  ax = fig.add_subplot(111)

  # Hide axes
  ax.axis('off')
  
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
        
  the_table = ax.table(cellText=clust_data,colLabels=collabel,loc='center')

  fig.canvas.set_window_title("Model num " + str(model))
  #for some reason zorder is not a keyword in ax.table
  the_table.set_zorder(10)

  global interactive
  if interactive:
    plt.show()
  
  if (not filename is None) :
    plt.savefig(filename + ".txt")


'''
Acceptable margin
To plot the acceptable intervalle of firing rates and the obtained FR
'''
def plot_acceptable_margin () :
  def labeling(xy, h, text):
    y = xy[1] + h + 1  # shift y-value for label so that it's above the rect
    plt.text(xy[0] + 0.05, y, text, ha="center", family='sans-serif', size=14)
  
  global NUCLEI,FRRNormal
  
  rect_size = 0.1
  nbNuclei = len(NUCLEI)
  
  fig=plt.figure(figsize=(rect_size*nbNuclei*10, 10))
  ax = fig.add_subplot(111)#, aspect='equal')
  
  # getting the ordinate max range
  ymax = 10
  
  # to boxplots
  #margin_data = list()
  
  for i,N in enumerate(NUCLEI) :
    #margin_data.append([FRRNormal[N][0],FRRNormal[N][1]])
  
    # saving the y range
    if (ymax < FRRNormal[N][1]) :
      ymax = FRRNormal[N][1] + 10 # with some margin
      
    x = 2*i*rect_size+rect_size
    y = FRRNormal[N][0]
    w = rect_size
    h = FRRNormal[N][1]-FRRNormal[N][0]
    
    # drawing a rectangle as the acceptable intervalle
    ax.add_patch(
        patches.Rectangle(
            # letting some margin
            (x, y),           # (x,y) position of the bottom left
            w,                # width
            h,                # height
            fill=False,       # remove background
        )
    )
    # setting a label on the rectangle
    labeling([x,y],h,N)
    
    # for each generated point for this nucleus
    # displaying and checking whether it is in or out side of the margins
    # TODO : getting back the generated points
    
    # To change to boxplots    
    #ax2.boxplot(margin_data)

  # removing labels from x
  ax.set_xticklabels([])

  ax.set_ylim([0,ymax])
  fig.canvas.set_window_title("Firing Rates margin")
  
  plt.show()





### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''
#plot_acceptable_margin()