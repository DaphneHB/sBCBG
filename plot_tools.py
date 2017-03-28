# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:32:40 2017

@author: daphnehb
"""
import os
import re
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from matplotlib import colors as mcolors

from LGneurons import NUCLEI, nbSim, FRRNormal, dataPath, FRRAnt, recType


### FUNCTIONS
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
Verifying if the current lineFR can be plotted
Meaning if it use the antagonist in antag
Return None or (antag string, [FRs])
'''
def can_plot(lineFR, antag, norm) :
  # the antagonist specifications
  antN, antInj = map(lambda x : re.sub(' ','',x), lineFR[1:3])
  antStr = antN + "_" + antInj.rstrip()
  if (antag is None or antag == 'all') and ('none' in antN) :
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
  
  # getting the ordinate max range
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
      plt.text(x + 0.05, -12, N, ha="center", family='sans-serif', size=10)
      
  else : # a specific antag
    # showing the GPi/e box only with a certain y margin
    antN, antInj = antag.split("_") # the antag string is Nucleus_injection form
    mnFR,mxFR = FRRAnt[antN][antInj]
    if (ymax < mxFR) :
      ymax = mxFR
    ## plotting every nuclei name
    for i,N in enumerate(NUCLEI) :
      x = 2*i*rect_size+rect_size
      if (N==antN) :
        # plotting the box for the injection site
        h = mxFR - mnFR
        # drawing a rectangle as the acceptable intervalle
        ax.add_patch(
            patches.Rectangle(
                # letting some margin
                (x, mxFR),           # (x,y) position of the bottom left
                rect_size,                # width
                h,                # height
                fill=False,       # remove background
            )
          )
        ymax = mxFR + h
      # setting a label on the rectangle
      plt.text(x + 0.05, -12, N, ha="center", family='sans-serif', size=14)
    ymax += 10
  # removing labels from x
  ax.set_xticklabels([])

  return ymax
  

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
    rnd = random.random()/10
    x_tab = np.arange(rect_size,xmax,rect_size*2) + rnd
    #print "X = ",x_tab," Y = ", listFR,"color = ",color
      
    ax.scatter(x_tab,listFR, c=color,marker=mark)
    # getting the maximum y possible :
    mx = max(listFR)
    if (ymax<mx) :
      ymax = mx
    
  legnd = []
  if model is None :
    mod_patches = []
    for m,cl in model_color.items() :
      mod_patches.append(patches.Patch(color=cl,label='model '+str(m)))
    legnd += mod_patches
  if not antag is None :
    ant_lgd = []
    for a,mk in anta_shape.items() :
      ant_lgd.append(mlines.Line2D([],[],color='black', marker=mk,markersize=10,label='ant ' + str(a) ))
    legnd += ant_lgd
  return legnd,ymax
  
'''
Works on the assumption that every simu firing rates
are reported in a allFiringRates.csv global file in the global log directory
the norm param decide whether or not we should have a normalization of the y axis
'''
def plot_margins_and_simus(filename=None,norm=False, antag=None, model=None) :
  global NUCLEI
  
  # antag & norm isnt possible
  if antag :
    norm = False  # setting norm to false auto
  
  # retrieving data in the input file allFiringRates.csv
  if filename is None :
    filename = 'allFiringRates'
  allFRfile = open('log/' + filename + ".csv",'r')
  allFRdata = allFRfile.readlines()
  allFRfile.close()
  # retrieving only the simu with the choosen model if there is
  if (not model is None) :
    allFRdata = filter(lambda x : ("#" + str(model)) in x ,allFRdata)
  if allFRdata == [] :
    print "---------- ERROR : No corresponding model simulated"
    return 1
  
  # getting the plot figure to fill
  fig = plt.figure()#figsize=(rect_size*nbNuclei*10, 5))
  ax = fig.add_subplot(111)
  
  nbNuclei = len(NUCLEI)

  rect_size = 0.1
  xmax = nbNuclei * 2 * rect_size        # at the end of the loop, x will get the last
  ymax = 0
  ax.set_xticklabels('')
  #### if we dont want a normalized plot lets draw the limit boxes
  if not norm :
    ymax = plot_margin_boxes(ax, rect_size, antag)
    ax.set_xticklabels([])
    ax.set_ylabel("Firing Rates (Hz)")
  # otherwise (normalized) lets label the x axis
  else :
    ymax = 1
    # Customize minor tick labels
    N_labels = ['' if e % 2 == 0 else NUCLEI[e/2] for e in range(nbNuclei * 2)] 
    ax.set_xticks(np.arange(0,xmax,rect_size), minor=True)
    ax.set_xticklabels(N_labels,minor=True)
    ax.set_yticks([0,1])
    ax.set_yticklabels(["Min","Max"])  
    ax.set_ylabel("Relative Margin")
  
  lgnd,ymax = plot_simu_points(allFRdata, ax, norm, antag, model, rect_size, [xmax,ymax])
  #### parametrizing the plot
  ax.legend(handles=lgnd,loc='upper center',bbox_to_anchor=(0.5,1.),ncol=3,fontsize='x-small').draggable()
  # legend x and y axis  
  ax.set_xlabel("BG Nuclei")
    
  ax.set_xlim([0,xmax + rect_size])       # with some margin  
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
  ax.grid()
  fig.canvas.set_window_title("Firing Rates margin" + title)
  
  # showing plot
  plt.show()

'''
Plotting number of simulations which get a score > score param
for each value of variable for the nucleus
nucleus is in NUCLEI
varible is either Ie or G
'''
def plot_score_ratio(variable, nucleus, dataPath=os.getcwd(), score=0, model=None) :
  NUM_COL = 0
  global NUCLEI
  if not nucleus in NUCLEI :
    print "------------ ERROR : Wrong nucleus"
    return 1
  if not ("Ie"==variable or "G"==variable) :
    print "------------ ERROR : Wrong variable name [" + variable + "]"
    return 1
   
  val_tab = []
  score_col = {}
  plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib
  n_var = variable + nucleus
  varN_values = {}   # dict {score : {val : nb}}
  model_pattern = re.compile("LG14modelID.*:\ *(\d+).")
  paramVal_pattern = re.compile(n_var + ".*:\ *(\d+).*")
  ################################## TODO review
  ### generate plot 
  for fName in os.listdir(dataPath) :
    dirPath = os.path.join(dataPath,fName)
    if os.path.isdir(dirPath) and fName.startswith("2017") :
      with open("score.txt","r") as scoreFile :
        obt_score = float(scoreFile.readline().rstrip())
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
        val = int(paramVal_pattern.findall(filter(lambda x : paramVal_pattern.search(x), Paramsdata)[0])[0])
        val_tab.append(val)
      except IndexError: # if there were no result : the variable name is wrong
        print "------------- ERROR : Wrong variable name [" + n_var + "]"
        return 1
      # extending the nb
      if (not varN_values.has_key(obt_score)) :
        varN_values[obt_score] = {val : 1} #dict(zip([float(x) for x in range(score,15)],[0.] * (15-score)))   # initializing every possible score for this value
      elif (varN_values[obt_score].has_key(val)) :
        varN_values[obt_score][val] += 1
      else :
        varN_values[obt_score][val] = 1
      # score color
      if (not score_col.has_key(obt_score)) :
        score_col[obt_score] = plot_colors[NUM_COL]
        
        
  # plot
  print varN_values
  width = 0.3
  plt_lgd = {}
  for obt_score,valsNb in varN_values.items() :
    # creating the unexisting vals ordering the list of vals
    val_list = []
    for v in val_tab :
      if not v in valsNb :
        val_list.append(0)
      else :
        val_list.append(valsNb[v])
    p = plt.bar(val_tab,val_list,width, color=score_col[obt_score])
    plt_lgd [p[0]] = obt_score
  # displaying
  plt.xlabel(n_var + " values")
  plt.ylabel('Number of simulations')
  plt.legend(plt_lgd.keys(),plt_lgd.values(), title="Scores")
  plt.title("Score of simulations according to the " + n_var + " param value")
  plot_margin = 0.25

  x0, x1, y0, y1 = plt.axis()
  plt.axis((x0 - plot_margin,
          x1 + plot_margin,
          y0 - plot_margin,
          y1 + plot_margin))
  plt.show()
  
  
  
### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''

plot_score_ratio("Ie","GPi",dataPath="data")