# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:32:40 2017

@author: daphnehb
"""
import os
import re
import random
import commands
import numpy as np
from operator import add

from pylab import cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from matplotlib import colors as mcolors

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.table import Table

import modelParams as mparams
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
        
  the_table = ax.table(cellText=clust_data,colLabels=collabel,loc='center',)

  fig.canvas.set_window_title("Model num " + str(model))
  plt.title("Model " + str(model))
  the_table.set_fontsize(20)
  #for some reason zorder is not a keyword in ax.table
  the_table.set_zorder(10)
  
  if (not filename is None) :
    plt.savefig(filename)
  else :
    plt.show()


'''
Plot the inDegree network graph for the given model
'''
def plot_inDegree_network(table_dict, model, filename=None) :
  pass



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
    NUM_COL = 0
    ant_col = {}
    plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib
    antNuclei = FRRAnt.keys()
    ## plotting every nuclei name
    for i,N in enumerate(NUCLEI) :
      x = 2*i*rect_size+rect_size
      # setting a label on the rectangle
      plt.text(x + 0.05, -12, N, ha="center", family='sans-serif', size=10)
      y = FRRNormal[N][0]
      h = FRRNormal[N][1]-FRRNormal[N][0]
      
      # drawing a rectangle as the acceptable intervalle
      ax.add_patch(
          patches.Rectangle(
              # letting some margin
              (x, y),           # (x,y) position of the bottom left
              rect_size,                # width
              h,                # height
              fill=False,       # remove background
          )
      )
      if N in antNuclei :
        # for every possible injection
        for antInj,boards in FRRAnt[N].items() :
          # getting the max possible value of y/x according to all antagonist injections
          mn,mx = boards
          # comparing max values
          if (ymax < mx) :
            ymax = mx
        ymax += 10
      
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
                (x, mnFR),           # (x,y) position of the bottom left
                rect_size,                # width
                h,                # height
                fill=False,       # remove background
            )
          )
        ymax = mxFR
      # setting a label on the rectangle
      plt.text(x + 0.05, -18, N, ha="center", family='sans-serif', size=10)
      
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
      
    ax.scatter(x_tab,listFR, c=color,marker=mark,edgecolor='')
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
def plot_margins_and_simus(filename=None,norm=False, antag=None, model=None, separated=None) :
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
Display on a plot the reason why a slope or bars cant ba displayed
'''
def plot_print_wrong(ax,reason) :
  if not ax is None :
    x0, x1, y0, y1 = ax.axis()
    ax.text(x1 - x0 + 2., y1 - y0, reason, ha="center", family='sans-serif', size=12)
  return 1
  
  
'''
Plotting number of simulations which get a score > score param
for each value of variable for the nucleus
nucleus is in NUCLEI
varible is either Ie or G
ATTENTION : plt.show() must be called after
'''
def plot_score_ratio(variable, nucleus, dataPath=os.getcwd(), score=0, model=None, axis=None,save=None) :
  NUM_COL = 0
  global NUCLEI
  if not nucleus in NUCLEI :
    reason = "------------ ERROR : Wrong nucleus"
    print reason
    return plot_print_wrong(axis,reason)
  if not ("Ie"==variable or "G"==variable) :
    reason = "------------ ERROR : Wrong variable name [" + variable + "]"
    print reason
    return plot_print_wrong(axis,reason)
  
  val_tab = []
  n_var = variable + nucleus
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
  # for every score
  for scKeys,valDict in varN_values.items() :
      # for values of the param that are not in score, put the number to 0
      for i,val in enumerate(val_tab) :
        if not valDict.has_key(val) :
          varN_values[scKeys][val] = 0
  
  # plot
  width = (max(val_tab)-min(val_tab))/(len(val_tab)*len(val_tab))
  plt_lgd = {}
  if axis is None :
    fig,ax = plt.subplots()
    fig.canvas.set_window_title("Score of simulations according to the " + n_var + " param value")
    ax.set_ylabel('Number of simulations')
  else :
    ax = axis
  
  # the bottom value of each hist
  # changing at each loop
  btm = [0] * len(val_tab)
  nbScores = len(varN_values)
  plot_colors = sorted(mcolors.cnames.values()[:nbScores + 1])     # list of colors in matplotlib
  for i,obt_score in enumerate(sorted(varN_values)) :
    # creating the unexisting vals ordering the list of vals
    valsNb = varN_values[obt_score]
    val_keys = sorted(valsNb)   # getting the key ordered
    val_nb = [valsNb[k] for k in val_keys]    # getting the value for the ordered keys
    p = ax.bar(val_keys,val_nb, width,color=plot_colors[i], bottom=btm,edgecolor='')
    plt_lgd [obt_score] = p[0]
    btm = map(add,btm,val_nb)
    
  mxVal = max(btm)    
  # displaying
  ax.set_yticks(np.linspace(0.,mxVal+1, 3 ) )
  ax.set_xlabel(n_var + " values")
  score_lgd = sorted(plt_lgd)
  ax.legend([plt_lgd[s] for s in score_lgd], score_lgd, title="Scores",loc=2,bbox_to_anchor=(0.9,1.1),fontsize='x-small').draggable()
  plot_margin = width

  x0, x1, y0, y1 = ax.axis()
  ax.axis((x0 - plot_margin,
          x1 + plot_margin,
          y0, 
          y1))
  ax.grid()
  if not save is None :
    fig.savefig(save)

'''
Plot the FR according to the variable G or Ie for a specific nucleus
val_tab is the list of tuple as (x,y,score)

for each score a specific color
'''
def plot_fr_by_var(n_var, val_tab, score_max, interv) :
  NUM_COLOR = 0
  score_col = {}
  plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib
  
  fig,ax = plt.subplots()
  # for (x,y,score)
  for prm, fr, sc in val_tab :
    if score_col.has_key(sc) :
      color = score_col[sc]
      ax.scatter([prm],[fr],c=color,edgecolor='')
      
    else :
      color = plot_colors[NUM_COLOR]
      score_col[sc] = color
      NUM_COLOR += 1
      ax.scatter([prm],[fr],c=color,label=sc,edgecolor='')
      
  ax.set_xticks(np.arange(*interv),minor=True)
  ax.set_xlabel(n_var + " values")
  ax.set_ylabel("Firing Rates")
  ax.legend(title="Scores (/" + str(score_max) + ")",loc='upper left',bbox_to_anchor=(0.,1.),ncol=3,fontsize='x-small').draggable()
  ax.set_title("Point cloud of " + str(n_var) + " with score")
  ax.grid()  
  
  fig.savefig("plots/FRby"+n_var+".png")
  plt.close(fig)
  
  
def plot_param_by_param(param1, param2, param3=None, dataPath=os.getcwd(), score=0, model=None, save=False) :
  NUM_COL = 0
  global NUCLEI
  # paramX_vals dict as : {param: (index in score_vals,)}
  param1_vals = []
  param2_vals = []
  param3_vals = []
  score_vals = []
  eachPoint = {}  # dict of list as point coordinate : score list
  
  fig = plt.figure(figsize=(8,6))
  if not param3 is None :
    axis = Axes3D(fig)
  else :
    axis = fig.add_subplot(111)
    
  model_pattern = re.compile("LG14modelID.*:\ *(\d+).")
  param1Val_pattern = re.compile(str(param1) + ".*:\ *(\d+[\.\d*]*).*")
  param2Val_pattern = re.compile(str(param2) + ".*:\ *(\d+[\.\d*]*).*")
  param3Val_pattern = re.compile(str(param3) + ".*:\ *(\d+[\.\d*]*).*")
  for fName in os.listdir(dataPath) :
    dirPath = os.path.join(dataPath,fName)
    if os.path.isdir(dirPath) and fName.startswith("2017") :
      with open(os.path.join(dirPath, "score.txt"),"r") as scoreFile :
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
  
  axis.set_xlabel(str(param1) + ' values')
  axis.set_ylabel(str(param2) + ' values')
  if param3 is None :
    # plotting 2D
    axis.scatter(param1_vals,param2_vals,c=cm.hsv(score_vals/max(score_vals)),s=500,marker='s',edgecolor='')
    title = "Score with x=" + str(param1) + " and y=" + str(param2)
    figname = param1 + "+" + param2 + "_score2D"
    axis.set_title(title)
    fig.canvas.set_window_title(figname)
    # displaying the score mean and median for each point
    for pt,scoreList in eachPoint.items() :
      meanScore = int(sum(scoreList) / (len(scoreList) * 14.) * 100 )
      medianScore = np.median(scoreList)
      txt = str(meanScore) + "%\n" +  str(medianScore)
      x,y = pt
      axis.text(x, y, txt, ha="center", family='sans-serif', size=8)
  else :
    axis.set_zlabel(str(param3) + ' values')
    # plotting 3D
    axis.scatter(param1_vals,param2_vals,param3_vals,c=cm.hsv(score_vals/max(score_vals)),s=500,marker='s',edgecolor='')
    title = "Score with x=" + str(param1) + ", y=" + str(param2) + " and z=" + str(param3)
    figname = param1 + "+" + param2 + "+" + param3 + "_score3D"
    axis.set_title(title)
    fig.canvas.set_window_title(figname)
    # displaying the score mean and median for each point
    for pt,scoreList in eachPoint.items() :
      meanScore = int(sum(scoreList) / (len(scoreList) * 14.) * 100 )
      medianScore = np.median(scoreList)
      txt = str(meanScore) + "%\n" +  str(medianScore)
      x,y,z = pt
      axis.text(x, y, z, txt, ha="center", family='sans-serif', size=8)
  cb = fig.colorbar(colmap)
  plt.subplots_adjust()
  if save:
    fig.savefig("log/" + figname + ".png")
  else :
    plt.show()
    
'''
Plotting for the 15 models for the 14 ranges, which one are/is wrong
for a given parametrization
'''
def plot_models_ranges(allFRdata, legend, paramFilePath=os.path.join(os.getcwd(),"modelParams.py"), models=np.arange(0,15,1),filename=None) :
  clust_data = []
  fig = plt.figure(1)
  
  ax = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=3)
  plt.title("Models'results for each range")
  ax.set_axis_off()  
  tb = Table(ax,bbox=[0,0,1,1])
  
  width = 0.2
  height = 1.0 / len(models)
  labels = []
  score = 0
  score_max = 0
  for row,model in enumerate(models) :
    frates = allFRdata[model]
    col = 0
    for frline in frates :
      frline = frline.split(',')[1:-1]
      # for each Nucleus, getting the results
      for nres in frline :
        score_max += 1
        nres = nres.rstrip().split('=')
        labels.append(nres[0])
        if nres[1]=="OK" :
          score += 1
          color = '#BFE8B7' # green
        else :
          color = '#E8B7B7' # red
        tb.add_cell(row,col,width,height,text=nres[1], loc='center',facecolor=color)
        col += 1
    tb.add_cell(row,col,width,height,text=str(score), loc='center',facecolor='white')
    col += 1
    break    
  labels.append("Score/" + str(score_max))
  # Row Labels...
  for i, label in enumerate(models):
      tb.add_cell(i, -1, width, height, text=label, loc='right', 
                  edgecolor='none', facecolor='none')
  # Column Labels...
  for j, label in enumerate(labels):
      tb.add_cell(len(models), j, width, height/5, text=label, loc='center', 
                         edgecolor='none', facecolor='none')
  ax.add_table(tb)
  
  # LEGEND
  rowSpan = max(3,1 / len(legend))
  ax = plt.subplot2grid((3,3), (0,2), colspan=1, rowspan=rowSpan)
  ax.set_axis_off()  
  plt.title("Parametrization used")
  tb = Table(ax,bbox=[0,0,1,1])
  height = 1. / len(legend)
  for row,lgd in enumerate(legend) :
    tb.add_cell(row,0,0.5,height,text=lgd,loc='left',edgecolor='white')
  ax.add_table(tb)
  
  fig.canvas.set_window_title("Passing LG14's tests")
  fig.tight_layout()
  fig.set_size_inches(w=11,h=7)
  if (not filename is None) :
    plt.savefig(filename)
  else :
    plt.show()

### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''
#plot_score_ratio("Ie","GPi",dataPath="/home/daphnehb/OIST/SangoTests/model2/copyBG")

#plot_models_ranges({0: ['#0 , MSN=OK , FSI=NO[+4.8340] , STN=OK , GPe=OK , GPi=NO[+6.6000] , ']},["'GMSN':5.7", "'GFSI':1.3", "'GSTN':1.38", "'GGPe':1.3", "'IeGPe':13.", "'GGPi':1.", "'IeGPi':11."],models=[0])