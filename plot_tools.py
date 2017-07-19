# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:32:40 2017

@author: daphnehb
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.table import Table
import scipy.interpolate as sp

from data_tools import *

### FUNCTIONS

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
  

def plot_param_legend(legend, plot_size,placement) : 
  # LEGEND
  rowSpan = max(plot_size[0]-1, 1 / len(legend))
  ax = plt.subplot2grid(plot_size, placement, colspan=plot_size[1] - placement[1], rowspan=rowSpan)
  ax.set_axis_off()  
  plt.title("Parametrization\nused")
  tb = Table(ax,bbox=[0,0,1,1])
  height = 1. / len(legend)
  for row,lgd in enumerate(legend) :
    tb.add_cell(row,0,1.,height,text=lgd,loc='left',edgecolor='white')
  ax.add_table(tb)
  
'''
To plot the table of the inDegree intervalle for a specific model
Must be called at the end of the simulation
'''
def plot_inDegrees_boarders_table(table_dict, model, filename=None) :
  # Labels for the column of the table
  collabel = ("Src", "Target", "Min inDegree", "Max inDegree", "Choosen Value")  
  clust_data = get_inDegree_to_plot(table_dict)
  
  nrows, ncols = len(clust_data)+1, len(collabel)
  hcell, wcell = 0.2, 1.
  hpad, wpad = 0, 1.5

  fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
  ax = fig.add_subplot(111)

  # Hide axes
  ax.axis('off')
  
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
allFiringRates is every lines satisfying model param
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
  
  allFRdata = get_data_from_file(lambda x : ("#" + str(model)) in x,filename=filename,model=model)
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
  n_var = variable + nucleus
  count = get_count_for_score(n_var, dataPath=dataPath, score=score, model=model,axis=axis)
  if type(count) is int :
    return count
  else :
    val_tab,varN_values = count
  # plot
  width = (max(val_tab)-min(val_tab))/(len(val_tab)*len(val_tab))
  plt_lgd = {}
  if axis is None :
    fig,ax = plt.subplots()
    fig.canvas.set_window_title("Score of simulations according to the " + n_var + " param value #"+str(model))
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
def plot_fr_by_var(n_var, val_tab, score_max, interv,model) :
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
  
  fig.savefig("plots/FRby"+n_var+"#"+str(model)+".png")
  plt.close(fig)
  
  
def plot_param_by_param(param1, param2, param3=None, dataPath=os.getcwd(), score=0, model=None, save=False) :
  NUM_COL = 0
  global NUCLEI
  # paramX_vals dict as : {param: (index in score_vals,)}
  
  fig = plt.figure(figsize=(8,6))
  if not param3 is None :
    axis = Axes3D(fig)
  else :
    axis = fig.add_subplot(111)
    
  param1_vals, param2_vals, param3_vals, eachPoint, score_vals, colmap = get_param_param_scores(param1, param2, param3=None, dataPath=dataPath, score=score, model=model)
  print "SCORES" ,score_vals
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
    figname = param1 + "+" + param2 + "+" + param3 + "_score3D_#" + str(model)
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
def plot_models_ranges(allFRdata, legend, models=np.arange(0,15,1),filename=None) :
  fig = plt.figure(1)
  plot_size = (8,10)
  
  ax = plt.subplot2grid(plot_size, (0,0), colspan=5, rowspan=8)
  plt.title("Models'results for each range")
  ax.set_axis_off()  
  tb = Table(ax,bbox=[0.1,0,1.7,1.])
  
  width = 0.25
  height = 1.0 / len(models)
  modLbl = []
  for modInd,model in enumerate(models) :
    frates = allFRdata[model]
    for frInd,frline in enumerate(frates) :
      row = modInd + frInd
      score = 0
      score_max = 0
      col = 0
      labels = ["Model"]
      modLbl.append(model)
      frline = frline.split(',')[1:-1]
      # for each Nucleus, getting the results
      for nres in frline :
        score_max += 1
        nres = nres.strip().split('=')
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
  labels.append("Score/" + str(score_max))
  
  # Row Labels...
  for i, label in enumerate(modLbl):
      tb.add_cell(i, -1, width, height, text=label, loc='right', 
                  edgecolor='none', facecolor='none')
  # Column Labels...
  global NUCLEI
  for j, label in enumerate(labels):
      if not (label in NUCLEI or label=="Model") :
        label = label.split("_")
        if len(label)==1 : # pour le score
          label = label[0].split("/")
          label[1] = "/" + label[1]
        label = label[0] + "\n" + label[1]
      tb.add_cell(len(modLbl), j-1, 0.2, height*2, text=label, loc='center', 
                         edgecolor='none', facecolor='none')
      '''ax.annotate(label,xy=(j*0.1,0),xycoords='axes fraction', ha='right',va='top',rotation=80,size=8)'''
  
  ax.add_table(tb)
  
  plot_param_legend(legend,plot_size,(1,8 ))
  
  fig.canvas.set_window_title("Passing LG14's tests")
  fig.tight_layout()
  fig.set_size_inches(w=11,h=7)
  if (not filename is None) :
    plt.savefig(filename)
  else :
    plt.show()


def plot_gap_from_range(vals, n_var, interv, nucleus_gap, model, param=None, filename=None) :
  plot_colors = mcolors.cnames.keys()[:len(nucleus_gap)]     # list of colors in matplotlib
  labels = []  
  
  if param is None:
    fig,ax = plt.subplots()
  else :
    fig = plt.figure()
    plot_size = (5,10)
    plot_param_legend(param,plot_size,(1,8))
    ax = plt.subplot2grid(plot_size, (0,0), colspan=7, rowspan=8)
  
  lower_gap = higher_gap = 0  
  # for each nucleus plotting its slope
  for ind,nucl in enumerate(nucleus_gap) :
    color = plot_colors[ind]
    higher_gap = max(higher_gap,max(nucleus_gap[nucl]))
    lower_gap = min(lower_gap,min(nucleus_gap[nucl]))
    labels.append(ax.plot(vals,nucleus_gap[nucl],'-^',c=color,label=color))
  ax.set_xticks(np.arange(*interv),minor=True)
  ax.set_xlabel(n_var + " values")
  upper_margin = abs(higher_gap - abs(lower_gap)) * 2
  ax.set_yticks(np.arange(lower_gap,higher_gap+upper_margin,1),minor=True)
  ax.set_ylabel("Firing Rates Gaps")
  ax.legend(nucleus_gap.keys(),title="Nuclei",loc='upper left',bbox_to_anchor=(0.,1.),ncol=3,fontsize='x-small').draggable()
  ax.set_title("FR gap slope for each nucleus with " + str(n_var) + " variations#" + str(model))
  ax.grid()

  fig.canvas.set_window_title("Firing Rates Gaps for model " +str(model))  
  if filename is None :
    plt.show()
  else :
    fig.savefig(filename)
    
    
def plot_score_by_value(parameter, pdata, simu_color, model=None, axis=None,filename=None) :
  # plot
  if axis is None :
    fig,ax = plt.subplots()
    fig.canvas.set_window_title("Best scores interpolation of simulations according to " + parameter + " #"+str(model))
    ax.set_ylabel('Score')
    ax.set_xlabel(parameter)
  else :
    ax = axis
  NUM_COLOR = 0
  plot_colors = mcolors.cnames.keys()     # list of colors in matplotlib  
  
  x = []
  y = []
  for val,simu_scoreList in pdata.items() :
    best_simu,best_score = simu_scoreList["best"]
    # getting the color for the simu
    if simu_color.has_key(best_simu) :
      color = simu_color[best_simu]
    else :
      color = plot_colors[NUM_COLOR]
      simu_color[best_simu] = color
      NUM_COLOR += 2
    # plotting
    ax.scatter([val],[int(best_score)], s=200, c=color,edgecolor='')
    x.append(val)
    y.append(int(best_score))
  # drawing a slope by interpolation
  '''
  interpolation = sp.interp1d(x,y,kind='cubic')
  ax.plot(xi,interpolation(xi))
  '''
  # plot adujstements
  ax.grid()
  xi = np.linspace(x[0],x[-1],(x[-1]-x[0])*10)
  ax.set_xticks(xi,minor=True)
  ax.set_yticks(y,minor=True)
  if not filename is None and axis is None :
    fig.savefig(filename)
  return simu_color
  
'''
Generate the plot with every piechart
Arguments :
  - xtab, the list of x values (also y value : symetric, square)
  - the data, a dict {value tuple : chan percentage tuple} where value tuple is the x,y coordinates
'''
def plot_multichan_pieChart(tab, values_dict, model,ratio,NbTrials,shuffled,nestseed,seed,save=None) :
  print "Plotting the 2-channels action selection competition for #%d" % model
  nbVals = float(len(tab)) + 1
  step = round(1/nbVals,1)
  cols=["grey","blue","green","red"]
  fig_range = [0.1,0.8]
  expe_range = [min(tab)-step,max(tab)+2*step]
  fig = plt.figure(figsize=(nbVals*0.9, nbVals*0.9))
  # main figure translated from 0.1 (0) to 0.8(1)
  ax1 = fig.add_axes([fig_range[0],fig_range[0],fig_range[1],fig_range[1]])
  rescale = lambda a : fig_range[0] + a * (fig_range[1] - fig_range[0]) + fig_range[0]*0.05   # assuming every experimental value is in between 0 and 1 + a small offset
  # specifications of the plot
  plt.xlabel('Channel 1\'s input activity')
  plt.ylabel('Channel 2\'s input activity')
  new_tab = np.arange(expe_range[0],expe_range[1],step)
  plt.xticks(new_tab)
  plt.yticks(new_tab)
  shuffleStr = "shuffled inputs" if shuffled else "non-shuffled inputs"
  plt.title('2-channels action selection competition\n#%d (over %d trials with %s)\n[ratio=%.3f, NEST seed=%d, random seed=%d]' % (model,NbTrials,shuffleStr,ratio,nestseed,seed))
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid()
  #legend
  p1 = Rectangle((0, 0), 1, 1, fc=cols[0])
  p2 = Rectangle((0, 0), 1, 1, fc=cols[1])
  p3 = Rectangle((0, 0), 1, 1, fc=cols[2])
  p4 = Rectangle((0, 0), 1, 1, fc=cols[3])
  plt.legend([p1,p2,p3,p4], ["No selection","Channel 1", "Channel 2", "Error"],bbox_to_anchor=(1.1, 1.05),fontsize='x-small')
  for x in tab :
    for y in tab :
      # aranging x, y for float equality comparison
      if not values_dict.has_key((x,y)):
        continue
      value_tuple = values_dict[(x,y)]
      # rescaling x and y to have them at the coordinates' intersection
      xpos = rescale(x)
      xpos = xpos + (1 - xpos) * 0.02
      ypos = rescale(y)
      ypos = ypos + (1 - ypos) * 0.02
      # inserting the piechart
      ax2 = fig.add_axes([xpos,ypos,step*0.7,step*0.7]) #fig.add_subplot(nbX,nbY,nb)
      ax2.axis("off")
      plot_piechart(ax2,value_tuple,colors=cols)
#  exit()
  # if the filename is not defined in save variable
  if save is None:
    print "\tShowing plot"
    plt.show()
  else :
    print "\tPlot saved under the name %s" % save
    fig.savefig(save)


'''
Generate the plot for a single piechart (for a single value tuple)
Parameters are :
  - the axes for this piechart (to insert it after with the others)
  - the percentages for everychannel (everytime in the same order)
  - the colors for each channel (everytime in the same order)
Return the axes of this plot (unused)
'''
def plot_piechart(axes,percentages,colors) :
  sizes = list(percentages)
  explode = [0] * len(percentages) 
  axes.pie(sizes, explode=explode,
          shadow=False, startangle=90,colors=colors) #, autopct='%2d%%')
  axes.set_aspect('equal', adjustable='box')
  return axes
  
def plot_fr_by_time(steps,firingRates,selected,actLevels,model,ratio,shuffled,nbTrials,simuTime,reversedChans,seeds,save=None) :
  print "Plotting the 2-channels action selection competition for #%d" % model
  cols=["grey","blue","green","red"]
  fig = plt.figure(figsize=(10,10))
  # specifications of the plot
  plt.xlabel('Time Steps',fontsize=12)
  plt.ylabel('Firing Rates (Hz)',fontsize=12)
  shuffleStr = "shuffled inputs" if shuffled else "non-shuffled inputs"
  revStr = "channels reversed" if reversedChans else ""
  seedStr = "seeds=[" + ",".join(map(str,seeds)) + "]"
  plt.title('2-channels action selection competition FR\n#%d (over %d trials with %s)\n[ratio=%.3f %s simuTime=%dms %s]\n\n' % (model,nbTrials,shuffleStr,ratio,seedStr,simuTime,revStr),fontsize=10)
  plt.grid()
  plt.xticks(steps,fontsize=10)
  ax2 = plt.twiny()
  actLevels = map(lambda x : "("+str(x[0])+",\n"+str(x[1])+")",actLevels)
  ax2.set_xticks(steps)
  ax2.set_xticklabels(actLevels,fontsize=7)
  #legend
  p1 = Rectangle((0, 0), 1, 1, fc=cols[0])
  p2 = Rectangle((0, 0), 1, 1, fc=cols[1])
  p3 = Rectangle((0, 0), 1, 1, fc=cols[2])
  p4 = Rectangle((0, 0), 1, 1, fc=cols[3])
  p5 = mlines.Line2D(range(1), range(1), color="white",marker='o',markersize=8,markerfacecolor="black")
  plt.legend([p1,p2,p3,p4,p5], ["No selection","Channel 1", "Channel 2", "Error","Selected"],bbox_to_anchor=(1.1, 1.),fontsize='x-small')
  # plotting
  avgFR1 = list()
  avgFR2 = list()
  for st in steps :
    thisStep = [[],[]]
    for sd in seeds :
      thisStep[0].append(firingRates[sd][0][st])
      thisStep[1].append(firingRates[sd][1][st])
    avgFR1.append(np.mean(thisStep[0]))
    avgFR2.append(np.mean(thisStep[1]))  
  plt.plot(steps,avgFR1,color=cols[1])
  plt.plot(steps,avgFR2,color=cols[2])
  for st,pt in zip(steps,selected) :
    plt.scatter([st],[0],s=500,c=cols[int(pt)],edgecolor='')

  # if the filename is not defined in save variable
  if save is None:
    print "\tShowing plot"
    plt.show()
  else :
    print "\tPlot saved under the name %s" % save
    fig.savefig(save)  

def draw_boxplot(ax,data, fill_color):
  bp = ax.boxplot(data, patch_artist=True,notch=False,vert=True)

  for patch in bp['boxes']:
      patch.set(facecolor=fill_color)       


def plot_errorFR_by_activity(actLevels,frtrials_dict,model,ratio,shuffled,nbTrials,simuTime,seed,reversedChans,save=None) :
  print "Plotting the 2-channels action selection competition for #%d" % model
  cols=["grey","blue","green","red"]
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_subplot(111)
  # specifications of the plot
  plt.xlabel('Input activities',fontsize=12)
  plt.ylabel('Firing Rates (Hz)',fontsize=12)
  #legend
  p2 = Rectangle((0, 0), 1, 1, fc=cols[1])
  p3 = Rectangle((0, 0), 1, 1, fc=cols[2])
  plt.legend([p2,p3], ["Channel 1", "Channel 2"],bbox_to_anchor=(1.1, 1.),fontsize='x-small')
  chan1vals = []
  chan2vals = []
  getValChan1 = lambda x : x[1]
  getValChan2 = lambda x : x[2]
  nbStep = 0
  # plotting the point one by one
  for key,dicVal in frtrials_dict.items() :
    # dicVal is a list of tuples (step,FR1,FR2,choosen channel)
    vals1 = map(getValChan1,dicVal)
    vals2 = map(getValChan2,dicVal)
    chan1vals.append(vals1)
    chan2vals.append(vals2)
    nbStep += 1
  # drawing channel 1
  draw_boxplot(ax,chan1vals,cols[1])
  # drawing channel 2
  draw_boxplot(ax,chan2vals,cols[2])

  actLevelStr = map(lambda x : "("+str(x[0])+",\n"+str(x[1])+")",actLevels)
  ax.set_xticklabels(actLevelStr,fontsize=7)
  plt.grid()
  shuffleStr = "shuffled inputs" if shuffled else "non-shuffled inputs"
  revStr = "channels reversed" if reversedChans else ""
  plt.title('2-channels action selection competition FR\n#%d means (%d simulations over %d trials with %s)\n[ratio=%.3f seed=%d simuTime=%dms %s]\n\n' % (model,nbStep,nbTrials,shuffleStr,ratio,seed,simuTime,revStr),fontsize=10)
  
  # if the filename is not defined in save variable
  if save is None:
    print "\tShowing plot"
    plt.show()
  else :
    print "\tPlot saved under the name %s" % save
    fig.savefig(save)
    
def plot_fr_for1(GPi_limits,oneChanFR,actLevels,model,offsetTime,simuTime,nestSeed,seeds,save=None) :
  print "Plotting the MC input activity test on 1 channel for #%d" % model
  fig = plt.figure(figsize=(10,10))
  # specifications of the plot
  plt.xlabel('Time Steps',fontsize=12)
  plt.ylabel('Firing Rates (Hz)',fontsize=12)
  seedStr = "seeds=NEST" + str(nestSeed) + "&[" + ",".join(map(str,seeds)) + "]"
  plt.title('1 channel FR analyze on MC plot\n#%d [%s offsetTime=%dms simuTime=%dms]\n\n' % (model,seedStr,offsetTime,simuTime),fontsize=10)
  plt.grid()
  actLevels = map(str,actLevels)
  steps = range(len(oneChanFR[seeds[0]]))
  plt.xticks(steps,fontsize=10)
  #legend
  p2 = Rectangle((0, 0), 1, 1, fc="blue")
  p4 = Rectangle((0, 0), 1, 1, fc="red")
  plt.legend([p2,p4], ["Channel 1","Gpi FR limits"],bbox_to_anchor=(1.1, 1.),fontsize='x-small')
  # plotting
  avgFR1 = list()
  for st in steps :
    thisStep = []
    for sd in seeds :
      thisStep.append(oneChanFR[sd][st])
    avgFR1.append(np.mean(thisStep))
  plt.plot(steps,avgFR1,color="blue")
  plt.plot(steps,[GPi_limits[0]] * len(steps),color='red')
  plt.plot(steps,[GPi_limits[1]] * len(steps),color='red')
  # if the filename is not defined in save variable
  if save is None:
    print "\tShowing plot"
    plt.show()
  else :
    print "\tPlot saved under the name %s" % save
    fig.savefig(save)
    
def plot_connectMap(nameSrc,nameTgt,nbTgtNeurons,nbChannels,specConnectMap,nPopulation,connectType,nestSeed,rndSeed,model=0,save=None) :
  plot_colors = mcolors.cnames.keys()[:nbChannels]     # list of colors in matplotlib
  #markers = [(2+i/2, 1+i%2, 0) for i in range(len(recType)*len(recType)*len(NUCLEI))]   # markers for the plot for 
  fig = plt.figure()
  ax = fig.add_subplot(111)
  xlbls = []
  for chan in range(nbChannels) :
    for neur in range(nbTgtNeurons) :      
      xVals = [nPopulation[chan][neur]]*len(specConnectMap[chan][neur])
      xlbls.append(int(xVals[0]))
      plt.scatter(xVals,list(specConnectMap[chan][neur]),c=plot_colors[chan])
  xlbls = [""] + sorted(xlbls)
  ax.set_xticklabels(xlbls)
  #plot_title = "Connection matrix of " + nameSrc + "->" + nameTgt + " in a " + str(nbChannels) + "-channel(s) simulation\nWith " + connectType + " connections #" + str(model)
  plot_title = "Connection matrix of " + nameSrc + "->" + nameTgt + " in a " + str(nbChannels) + "-channel(s) simulation\n#" + str(model) + " NEST seed=" + str(nestSeed) + ", random seed=" + str(rndSeed)
  plt.title(plot_title)
  print "Connect Map : ",specConnectMap
  # if the filename is not defined in save variable
  if save is None:
    print "\tShowing plot"
    plt.show()
  else :
    print "\tPlot saved under the name %s" % save
    fig.savefig(save)
    
  
#plot_errorFR_by_activity([(0.0, 0.0), (0.0, 0.5), (0.0, 0.0), (0.0, 0.5), (0.0, 0.0), (0.0, 0.5), (0.0, 0.0), (0.0, 0.5), (0.0, 0.0), (0.0, 0.5)],{(0.0, 0.0): [(0, 70.0, 65.238095238095241, '0'), (2, 73.571428571428569, 67.857142857142847, '0'), (4, 65.714285714285708, 63.095238095238095, '0'), (6, 63.095238095238095, 70.238095238095241, '0'), (8, 70.0, 68.095238095238102, '0')], (0.0, 0.5): [(1, 90.0, 10.0, '2'), (3, 89.761904761904759, 6.9047619047619051, '2'), (5, 84.285714285714278, 13.80952380952381, '2'), (7, 86.904761904761898, 11.428571428571429, '2'), (9, 93.095238095238088, 5.0, '2')]},9,1.5,False,5,300,False,None)
#plot_multichan_pieChart(np.arange(0,1.1,0.1),{(0.1,0.7):(0,10,90),(0.,0):(0,50,50),(0.6,0.5):(10,60,30)})


### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''
#plot_score_ratio("Ie","GPi",dataPath="/home/daphnehb/OIST/SangoTests/model2/copyBG")

#plot_models_ranges({0: ['#0 , MSN=OK , FSI=NO[+4.8340] , STN=OK , GPe=OK , GPi=NO[+6.6000] , ']},["'GMSN':5.7", "'GFSI':1.3", "'GSTN':1.38", "'GGPe':1.3", "'IeGPe':13.", "'GGPi':1.", "'IeGPi':11."],models=[0])

#plot_fr_by_time([0, 1, 2, 3, 4],{17: [[ 78.57142857,  65.23809524,  69.52380952,  64.28571429,73.33333333],[ 72.38095238,  64.76190476,  63.80952381,  60.,59.04761905]]},['0', '0', '0', '0', '3'],[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)],9,1.5,False,1,150,False,[17])

#plot_fr_for1([59.1,79.5],{17:{0: 70.15714285714286, 1: 70.47142857142858, 2: 70.77142857142857, 3: 70.02857142857142, 4: 71.34285714285714, 5: 68.74285714285715, 6: 70.74285714285715, 7: 70.61428571428571, 8: 72.44285714285714, 9: 69.91428571428573, 10: 68.88571428571429, 11: 69.27142857142857, 12: 71.6, 13: 69.81428571428572, 14: 69.32857142857144, 15: 69.07142857142857, 16: 70.44285714285714, 17: 69.45714285714286, 18: 70.55714285714286, 19: 70.0, 20: 69.75714285714285, 21: 70.92857142857143, 22: 70.57142857142858, 23: 68.9, 24: 69.62857142857143, 25: 70.12857142857143, 26: 69.34285714285714, 27: 69.62857142857143, 28: 69.2, 29: 71.12857142857143, 30: 72.11428571428571, 31: 70.64285714285715, 32: 70.27142857142857, 33: 71.21428571428572, 34: 69.61428571428571, 35: 69.54285714285714, 36: 69.07142857142857, 37: 71.75714285714285, 38: 70.57142857142858, 39: 69.4857142857143, 40: 69.11428571428571, 41: 70.1, 42: 70.34285714285714, 43: 69.58571428571429, 44: 69.81428571428572, 45: 68.12857142857143, 46: 70.38571428571429, 47: 70.74285714285715, 48: 71.08571428571427, 49: 69.44285714285714},1:{0: 71.24285714285715, 1: 70.4857142857143, 2: 69.7, 3: 70.52857142857142, 4: 70.15714285714286, 5: 69.82857142857144, 6: 69.9857142857143, 7: 70.3, 8: 70.7, 9: 68.8, 10: 68.44285714285714, 11: 70.25714285714285, 12: 69.64285714285714, 13: 70.78571428571429, 14: 70.54285714285714, 15: 69.8, 16: 70.28571428571428, 17: 68.24285714285715, 18: 70.2, 19: 71.18571428571428, 20: 68.94285714285714, 21: 70.07142857142857, 22: 70.88571428571429, 23: 71.11428571428571, 24: 71.25714285714285, 25: 69.85714285714285, 26: 68.72857142857143, 27: 68.84285714285714, 28: 71.21428571428572, 29: 69.68571428571428, 30: 70.14285714285714, 31: 68.55714285714286, 32: 69.9857142857143, 33: 71.68571428571428, 34: 69.21428571428571, 35: 70.0142857142857, 36: 69.52857142857142, 37: 72.44285714285714, 38: 70.17142857142856, 39: 69.4857142857143, 40: 69.8, 41: 70.81428571428572, 42: 69.15714285714286, 43: 69.77142857142857, 44: 69.38571428571429, 45: 68.15714285714286, 46: 71.3, 47: 70.41428571428573, 48: 71.64285714285715, 49: 70.6}},np.zeros(50),9,1000,5000,[1])

#plot_fr_for1([59.1, 79.5],{1: ((245, 5), {0: 71.24285714285715, 1: 70.4857142857143, 2: 69.7, 3: 70.52857142857142, 4: 70.15714285714286, 5: 69.82857142857144, 6: 69.9857142857143, 7: 70.3, 8: 70.7, 9: 68.8, 10: 68.44285714285714, 11: 70.25714285714285, 12: 69.64285714285714, 13: 70.78571428571429, 14: 70.54285714285714, 15: 69.8, 16: 70.28571428571428, 17: 68.24285714285715, 18: 70.2, 19: 71.18571428571428, 20: 68.94285714285714, 21: 70.07142857142857, 22: 70.88571428571429, 23: 71.11428571428571, 24: 71.25714285714285, 25: 69.85714285714285, 26: 68.72857142857143, 27: 68.84285714285714, 28: 71.21428571428572, 29: 69.68571428571428, 30: 70.14285714285714, 31: 68.55714285714286, 32: 69.9857142857143, 33: 71.68571428571428, 34: 69.21428571428571, 35: 70.0142857142857, 36: 69.52857142857142, 37: 72.44285714285714, 38: 70.17142857142856, 39: 69.4857142857143, 40: 69.8, 41: 70.81428571428572, 42: 69.15714285714286, 43: 69.77142857142857, 44: 69.38571428571429, 45: 68.15714285714286, 46: 71.3, 47: 70.41428571428573, 48: 71.64285714285715, 49: 70.6}), 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None},[ 0.  ,0. , 0. , 0.  ,0.  ,0.  ,0.  ,0. , 0. , 0. , 0. , 0. , 0.,  0.,  0.,  0. , 0.,  0.,  0.  ,0.  ,0.  ,0.  ,0. , 0. , 0.,  0.,  0. , 0.,  0.,  0.,  0.,  0.,  0.,  0. , 0. , 0.,  0.,  0. , 0.,  0.,  0.,  0. , 0.,  0. , 0. , 0.  ,0.  ,0. , 0. , 0.],9,1000,5000,[1],None)

#plot_connectMap("CMPf","STN",8,1,[[(8945, 8951, 8952, 8946, 8950, 8949, 8953, 8948, 8947), (8953, 8945, 8946, 8952, 8947, 8950, 8949, 8948, 8951), (8948, 8953, 8950, 8949, 8947, 8951, 8952, 8945, 8946), (8948, 8947, 8946, 8950, 8949, 8953, 8951, 8945, 8952), (8953, 8950, 8945, 8946, 8947, 8948, 8952, 8951, 8949), (8953, 8948, 8946, 8951, 8947, 8950, 8952, 8949, 8945), (8947, 8948, 8953, 8945, 8946, 8951, 8952, 8949, 8950), (8949, 8951, 8948, 8947, 8953, 8950, 8945, 8946, 8952)]],[[2698,2699,2700,2701,2702,2703,2704,2705]],"focused",model=9)