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
def plot_multichan_pieChart(tab, values_dict) :
  nbVals = float(len(tab)) + 1
  step = round(1/nbVals,1)
  cols=["grey","yellow","red"]
  fig_range = [0.1,0.8]
  expe_range = [min(tab)-step,max(tab)+2*step]
  fig = plt.figure(figsize=(nbVals*0.9, nbVals*0.9))
  # main figure translated from 0.1 (0) to 0.8(1)
  ax1 = fig.add_axes([fig_range[0],fig_range[0],fig_range[1],fig_range[1]])
  rescale = lambda a : fig_range[0] + a * (fig_range[1] - fig_range[0]) + fig_range[0]*0.05   # assuming every experimental value is in between 0 and 1 + a small offset
  # specifications of the plot
  plt.xlabel('Channel 1')
  plt.ylabel('Channel 2')
  new_tab = np.arange(expe_range[0],expe_range[1],step)
  plt.xticks(new_tab)
  plt.yticks(new_tab)  
  plt.title('2-channels action selection competition')
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid()
  #legend
  p1 = Rectangle((0, 0), 1, 1, fc=cols[0])
  p2 = Rectangle((0, 0), 1, 1, fc=cols[1])
  p3 = Rectangle((0, 0), 1, 1, fc=cols[2])
  plt.legend([p1,p2,p3], ["No selection","Channel 1", "Channel 2"],bbox_to_anchor=(1.1, 1.05),fontsize='x-small')
  for x in tab :
    for y in tab :
      # aranging x, y for float equality comparison
      x = round(x,1)
      y = round(y,1)
      if not values_dict.has_key((x,y)):
        continue
      value_tuple = values_dict[(x,y)]
      # rescaling x and y to have them at the coordinates' intersection
      x = rescale(x)
      x = x + (1 - x) * 0.03
      y = rescale(y)
      y = y + (1 - y) * 0.03
      # inserting the piechart
      ax2 = fig.add_axes([x,y,step*0.7,step*0.7]) #fig.add_subplot(nbX,nbY,nb)
      ax2.axis("off")
      plot_piechart(ax2,value_tuple,colors=cols)
#  exit()
  plt.show()


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

#plot_multichan_pieChart(np.arange(0,1.1,0.1),{(0.1,0.7):(0,10,90),(0.,0):(0,50,50),(0.6,0.5):(10,60,30)})


### Tests
'''
table = {'MSN->GPe': (105.37051792828686, 18018.358565737053), 'MSN->GPi': (151.65986013986014, 31696.91076923077), 'GPe->GPi': (1.4744055944055943, 23.59048951048951), 'GPe->MSN': (0.0015184513006654568, 0.14121597096188748), 'GPe->GPe': (0.84, 31.919999999999998), 'CMPf->GPe': (0.3426294820717131, 15.760956175298803), 'CMPf->GPi': (0.6013986013986014, 83.59440559440559), 'CMPf->FSI': (0.16165413533834586, 122.21052631578947), 'PTN->FSI': (-1, 5.0), 'CMPf->STN': (1.1168831168831168, 64.77922077922078), 'STN->MSN': (0.0004949334543254689, 0.05394774652147611), 'GPe->STN': (3.25974025974026, 61.935064935064936), 'STN->GPe': (0.2546215139442231, 74.34948207171315), 'STN->GPi': (0.38769230769230767, 63.96923076923076), 'CMPf->MSN': (0.003251663641863279, 7.244706594071385), 'FSI->FSI': (1.0, 140.0), 'CSN->MSN': (-1, 318.0), 'PTN->MSN': (-1, 5.0), 'FSI->MSN': (0.020114942528735632, 43.689655172413794), 'MSN->MSN': (1.0, 509.0), 'PTN->STN': (-1, 262.0), 'CSN->FSI': (-1, 489.0), 'GPe->FSI': (0.07548872180451129, 36.46105263157895)}

plot_inDegrees_boarders_table(table,'0')
'''
#plot_score_ratio("Ie","GPi",dataPath="/home/daphnehb/OIST/SangoTests/model2/copyBG")

#plot_models_ranges({0: ['#0 , MSN=OK , FSI=NO[+4.8340] , STN=OK , GPe=OK , GPi=NO[+6.6000] , ']},["'GMSN':5.7", "'GFSI':1.3", "'GSTN':1.38", "'GGPe':1.3", "'IeGPe':13.", "'GGPi':1.", "'IeGPi':11."],models=[0])

