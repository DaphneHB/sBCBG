# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:32:40 2017

@author: daphnehb
"""
import re
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors as mcolors

from LGneurons import NUCLEI, nbSim, FRRNormal, dataPath, FRRAnt, recType


### FUNCTIONS
'''
According to the norm param, the function will normalize 
from certain ranges in FRRNormal to 0-1 the input in frlist.
'''
def normalize(frlist,norm=False) :
  if (not norm or len(n_boarders)!=len(frlist)) :
    return frlist
  new_list = [0] * 5
  print "normalizing"
  global NUCLEI
  for ind,N in enumerate(NUCLEI) :
    new_list[ind] = (frlist[ind]-FRRNormal[N][0])/float(FRRNormal[N][1] - FRRNormal[N][0])
  return new_list
  
'''
Verifying if the current lineFR can be plotted
Meaning if it use the antagonist in antag
Return None or (antag string, [FRs])
'''
def can_plot(lineFR, antag, norm) :
  # the antagonist specifications
  antN, antInj = lineFR[1:3]
  antStr = antN + "_" + antInj

  if (antag is None and 'none' in antN) :
    lineFR = normalize(map(float,lineFR[3:-1]),norm)  # firing rates values
    antStr = 'none'
  elif (antag == 'all' and not 'none' in antN) or (antag == antStr):
    lineFR = normalize(map(float,lineFR[3:-1]),norm)  # firing rates values
  else :
    return None
    
  return antStr, lineFR
  

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

  plt.show()
  
  if (not filename is None) :
    plt.savefig(filename + ".txt")


'''
Acceptable margin for a specific model
To plot the acceptable intervalle of firing rates and the obtained FR
'''
def plot_acceptable_margin_oneModel (model=None) :
  
  def labeling(xy, h, text):
    y = xy[1] + h + 1  # shift y-value for label so that it's above the rect
    plt.text(xy[0] + 0.05, y, text, ha="center", family='sans-serif', size=14)
    
  def column(matrix, i):
    return [float(row[i]) for row in matrix]    # from string list to float list
  
  global NUCLEI,FRRNormal
  
  simusFR = list()
  # getitng the different existing firing rates in the log/firingRates.csv file
  firingRatesFile=open(dataPath+'firingRates.csv','r')
  allFiringRates = firingRatesFile.readlines()
  firingRatesFile.close()
  # for each simulation done for this specific model
  # recording in simusFR table the firing rates
  for lineSimu in allFiringRates :
    simuCases = lineSimu.split(",")
    # if it is the good model
    if (("#" + str(model)) in simuCases[0]) :
      # for no antagonnist injection (normal)
      if ("none" in simuCases[1]) :
        simusFR.append(simuCases[2:-1])     # removing the last \n
  nbSimus = len(simusFR)
    
  rect_size = 0.1
  nbNuclei = len(NUCLEI)
  
  fig=plt.figure(figsize=(rect_size*nbNuclei*10, 10))
  ax = fig.add_subplot(111)
  
  # getting the ordinate max range
  ymax = -10
  
  # to boxplots
  #margin_data = list()
  
  for i,N in enumerate(NUCLEI) :
    #margin_data.append([FRRNormal[N][0],FRRNormal[N][1]])
  
    # saving the y range
    if (ymax < FRRNormal[N][1]) :
      ymax = FRRNormal[N][1]
      
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
    
    # for each simulation done and registered
    # plotting the points
    # getting some random x around a certain value (for legibility)
    xPt = np.random.normal(x + w/2, rect_size/5, nbSimus)       # mean loc, scale/variance, nb
    # getting the corresponding column
    yPt = column(simusFR,i)
    ax.scatter(xPt,yPt,s=10)
    # to verify the maximum value to plot even with the points
    if (max(yPt) > ymax) :
      ymax = max(yPt)
    
    # To change to boxplots    
    #ax.boxplot(margin_data)

  # removing labels from x
  ax.set_xticklabels([])

  ax.set_ylim([-5,ymax + 10])       # with some margin
  fig.canvas.set_window_title("Firing Rates margin for model " +str(model))
  
  plt.show()

'''
Add the plotted margin boxis to the axis ax according to antag inj
If there is no antag inj : every boxes are plot
if we want to plot a specific antag inj only one box is plotted
if we want to plot every antag inj : no box is plot
Return the global xmax and ymax
'''
def plot_margin_boxes(ax, rect_size, antag) :
  
  def labeling(xy, h, text):
    y = xy[1] + h + 1  # shift y-value for label so that it's above the rect
    plt.text(xy[0] + 0.05, y, text, ha="center", family='sans-serif', size=14)
    
  global NUCLEI,FRRNormal, FRRAnt
  
  rect_size = 0.1
  nbNuclei = len(NUCLEI)

  # getting the ordinate max range
  xmax = 0
  ymax = 0
  
  # to boxplots
  #margin_data = list()
  if (antag is None) :
    ## plotting every box margin
    for i,N in enumerate(NUCLEI) :
      #margin_data.append([FRRNormal[N][0],FRRNormal[N][1]])
    
      # saving the y range
      if (ymax < FRRNormal[N][1]) :
        ymax = FRRNormal[N][1]
        
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
      
      # To change to boxplots
      #ax.boxplot(margin_data)
  elif (antag == "all") :
    # getting the max possible value of y/x according to all antagonist injections
    for antN in FRRAnt.keys() :
      # for every possible injection
      for mn,mx in FRRAnt[antN].values() :
        # comparing max values
        if (ymax < mx) :
          ymax = mx
    ymax += 10
    ## plotting every nuclei name
    for i,N in enumerate(NUCLEI) :
      x = 2*i*rect_size+rect_size      
      # setting a label on the rectangle
      plt.text(x + 0.05, -8, N, ha="center", family='sans-serif', size=10)
      
  else : # a specific antag
    # showing the GPi/e box only with a certain y margin
    antN, antInj = antag.split("_") # the antag string is Nucleus_injection form
    mn,mx = FRRAnt[antN][antInj]
    if (ymax < mx) :
      ymax = mx
    ## plotting every nuclei name
    for i,N in enumerate(NUCLEI) :
      x = 2*i*rect_size+rect_size
      if (N==antN) :
        # plotting the box for the injection site
        h = mx - mn
        # drawing a rectangle as the acceptable intervalle
        ax.add_patch(
            patches.Rectangle(
                # letting some margin
                (x, mx),           # (x,y) position of the bottom left
                rect_size,                # width
                h,                # height
                fill=False,       # remove background
            )
          )
      # setting a label on the rectangle
      plt.text(x + 0.05, -10, N, ha="center", family='sans-serif', size=14)
    ymax += 10
  xmax = nbNuclei * 2 * rect_size        # at the end of the loop, x will get the last
  # removing labels from x
  ax.set_xticklabels([])

  return (xmax,ymax)

'''
Acceptable margin for every tested model which appear in firingRates.csv
To plot the acceptable intervalle of firing rates and the obtained FR
'''
def plot_acceptable_margin_all () :
  NUM_COLOR = 1
  
  def labeling(xy, h, text):
    y = xy[1] + h + 1  # shift y-value for label so that it's above the rect
    plt.text(xy[0] + 0.05, y, text, ha="center", family='sans-serif', size=14)
    
  def column(matrix, i):
    return [float(row[i]) for row in matrix]    # from string list to float list
  
  global NUCLEI,FRRNormal
  
  rect_size = 0.1
  nbNuclei = len(NUCLEI)
  
  fig=plt.figure()#figsize=(rect_size*nbNuclei*10, 5))
  ax = fig.add_subplot(111)
  
  # getting the ordinate max range
  xmax = 0
  ymax = 0
  
  # to boxplots
  #margin_data = list()
  
  ## plotting every box margin
  for i,N in enumerate(NUCLEI) :
    #margin_data.append([FRRNormal[N][0],FRRNormal[N][1]])
  
    # saving the y range
    if (ymax < FRRNormal[N][1]) :
      ymax = FRRNormal[N][1]
      
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
    xmax = x+rect_size        # at the end of the loop, x will get the last
    
    # To change to boxplots    
    #ax.boxplot(margin_data)

  # getitng the different existing firing rates in the log/firingRates.csv file
  firingRatesFile=open(dataPath+'firingRates.csv','r')
  allFiringRates = firingRatesFile.readlines()
  firingRatesFile.close()
  
  model_color = {}              # dico[model] = color for plot
  plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib
  models_labels = []    # to remember which model was labeled
  
  ## plotting every simu (generating a color by simu)
  for lineSimu in allFiringRates :
    lineSimu = lineSimu.split(",")       # from a line string to a string array
    simu_num = int(re.findall('\d+',lineSimu[0])[0])
    print "Simu num : ",simu_num
    # getting the (new) color for this model
    if (model_color.has_key(simu_num)) :
        color = model_color[simu_num]
    else :
        color = plot_colors[NUM_COLOR]
        model_color[simu_num] = color
        NUM_COLOR += 1
    # retreiving the normal simu
    if ('none' in lineSimu[1]) :
      lineSimu = map(float,lineSimu[2:-1])
      rnd = random.random() + 1
      x_tab = np.linspace(rect_size*rnd, xmax/rnd, nbNuclei)
      print "X = ",x_tab," Y = ", lineSimu,"color = ",color
      if not simu_num in models_labels :
        models_labels.append(simu_num)
        ax.scatter(x_tab,lineSimu, c=color,label='model '+str(models_labels),marker='s')
      else :
        ax.scatter(x_tab,lineSimu, c=color,marker='s')
    # getting the maximum y possible :
    mx = max(lineSimu)
    if (ymax<mx) :
      ymax = mx
      
  # removing labels from x
  ax.set_xticklabels([])
  
  ax.legend()
  
  ax.set_xlim([0,xmax + rect_size*10])       # with some margin  
  ax.set_ylim([-5,ymax + 10])       # with some margin
  fig.canvas.set_window_title("Firing Rates margin for multiple models ")
  
  plt.show()

'''
Plot every point of each simulation satisfying the antag and model params
allFiringRates is every lines satisfying model
ax the plot axis
norm : whether the plot is normalized or not
antag define the antagonist results showing
model is the model number
rect_size is the size or a nucleus representation in the plot
xyMax is a list (xmax,ymax)
Return ymax value
'''
def plot_simu_points(allFiringRates, ax, norm, antag, model, rect_size, xyMax) :
  NUM_COLOR = 0
  NUM_MARK = 1

  def column(matrix, i):
    return [float(row[i]) for row in matrix]    # from string list to float list
  
  global NUCLEI, FRRNormal

  xmax,ymax = xyMax
  nbNuclei = len(NUCLEI)

  model_color = {}              # dico[model] = color for plot
  # if we want to plot every models
  if (model is None) : 
    plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib
  else :
    model_color[int(model)] = color = 'blue'
  models_labels = []    # to remember which model was labeled
  if (antag == 'all') :
    markers = [(2+i/2, 1+i%2, 0) for i in range(len(recType)*len(recType)*len(NUCLEI))]   # markers for the plot for antag
    anta_shape = {'none' : markers[0]}   # initializing a shape for no ant inj
  elif (antag is None) :
    anta_shape = {'none' : 'o'}
    # by default
    mark = anta_shape['none']
  else :
    anta_shape = {antag : '^'}
    # by default
    mark = anta_shape[antag]


  for lineSimu in allFiringRates :
    lineSimu = lineSimu.split(",")       # from a line string to a string array

    # either None or (antag, list FR)
    can_plot_res = can_plot(lineSimu,antag,norm)
    if (can_plot_res is None) :
      continue
    else :
      antStr, listFR = can_plot_res
    
    # the model number
    model_num = int(re.findall('\d+',lineSimu[0])[0])
    if (model is None) :
      # getting the (new) color for this model
      if (model_color.has_key(model_num)) :
        color = model_color[model_num]
      else :
        color = plot_colors[NUM_COLOR]
        model_color[model_num] = color
        NUM_COLOR += 1
    # register shapes for every antag
    if (antag == 'all') :
      if (anta_shape.has_key(antStr)) :
        mark = anta_shape[antStr]
      else :
        mark = markers[NUM_MARK]
        anta_shape[antStr] = mark
        NUM_MARK += 1
    
    # x list coordinates
    rnd = random.random()/10 + 1
    x_tab = np.arange(rect_size,xmax,rect_size*2) * rnd #np.linspace(rect_size*rnd, xmax/rnd, nbNuclei)
    print "X = ",x_tab," Y = ", listFR,"color = ",color
    
    if not model_num in models_labels :
      models_labels.append(model_num)
      ax.scatter(x_tab,listFR, c=color,label='model '+str(models_labels),marker=mark)
    else :
      ax.scatter(x_tab,listFR, c=color,marker=mark)
    # getting the maximum y possible :
    mx = max(lineSimu)
    if (ymax<mx) :
      ymax = mx
    print ymax, type(ymax)
  return ymax
  
'''
Works on the assumption that every simu firing rates
are reported in a allFiringRates.csv global file in the global log directory
the norm param decide whether or not we should have a normalization of the y axis
'''
def plot_margins_and_simus(norm=False, antag=None, model=None) :
  global NUCLEI
  
  # antag & norm isnt possible
  if antag :
    norm = False  # setting norm to false auto
  
  # getting the plot figure to fill
  fig = plt.figure()#figsize=(rect_size*nbNuclei*10, 5))
  ax = fig.add_subplot(111)
  
  rect_size = 0.1
  xmax = ymax = 0
  #### if we dont want a normalized plot lets draw the limit boxes
  if not norm :
    xmax, ymax = plot_margin_boxes(ax, rect_size, antag)
    ax.set_ylabel("Firing Rates (Hz)")
  # otherwise lets label the x axis
  else :
    ax.set_xticklabels(NUCLEI)
    ax.set_yticklabels(["Min","Max"])
    ax.set_ylabel("Relative Margin")

  #### plotting simu points
  # retrieving data in the input file allFiringRates.csv
  allFRfile = open('log/allFiringRates.csv','r')
  allFRdata = allFRfile.readlines()
  allFRfile.close()
  # retrieving only the simu with the choosen model if there is
  if (not model is None) :
    allFRdata = filter(lambda x : ("#" + str(model)) in x ,allFRdata)
  
  ymax = plot_simu_points(allFRdata, ax, norm, antag, model, rect_size, [xmax,ymax])
  #### parametrizing the plot
  ax.legend()
  # legend x and y axis  
  ax.set_xlabel("BG Nuclei")
    
  ax.set_xlim([0,xmax + rect_size*10])       # with some margin  
  if (not norm) :
    ax.set_ylim([-5,ymax + 10])       # with some margin
  # setting the name according to the params
  title = ""
  if (norm) :
    title += " normalized [0-1]"
  if (model is None) :
    title += " (*models)"
  else :
    title += " (model" + str(model) + ")"
  if (not antag is None) :
    title += " - antagonist injection"
  
  fig.canvas.set_window_title("Firing Rates margin" + title)
  
  # showing plot
  #plt.show()

### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''

#plot_acceptable_margin_all()
plot_margins_and_simus(norm=False, model=9)