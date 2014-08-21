# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to analyze and plot data series

This module defines the following public functions:
find_equilibration
parse_series
plot_series
write_series

Can be executed from the command line as a stand-alone program
"""

import os
import logging

import numpy as np
import matplotlib.pylab as plt
from scipy.stats import spearmanr
from scipy.stats import kendalltau

import simulationobjects

logger = logging.getLogger('protoms')

def find_equilibration(x,y,atleast=50,threshold=0.05,nperm=0) :
  """
  Find the equilibration time of a data series

  Parameters
  ----------
  x : Numpy array
    the x data
  y : Numpy array
    the y data
  atleast : int, optional
    this is the smallest number of snapshots considered to be production
  threshold : float, optional
    the confidence level
  nperm : int, optional
    if greater than zero, a permutation test is peformed with nperm synthetic permuations

  Returns
  -------
  int
    the length of the equilibration period
  """
  for i in range(0,y.shape[0]-atleast) :
    if nperm == 0 : # Performs an assymptotic test, fast
      tau,p = kendalltau(x[i:],y[i:])
      if p > threshold : return i
    else : # Performs a rigorous permutation test, slow
      rhos = np.zeros(nperm)
      rho0,p = spearmanr(x[i:],y[i:])
      for j in range(nperm) :
        rhos[j],p = spearmanr(x[i:],np.random.permutation(y[i:]))
      ncnt = (rho0*rho0 < rhos*rhos).sum()
      if float(ncnt) / float(nperm) > threshold : return i
  return y.shape[0]-atleast

def parse_series(series,results) :
  """
  Extract data series and label from a results file based on a label
  
  Parameters
  ----------
  series : list of string
    the series to extract
  results : SnapshotResults object
    all the results
  
  Returns
  -------
  list of NumpyArrays
    the data series
  list of string
    the labels for the series 
  """

  def parse_compound_series(series,results) :
    """
    Parse a compound series, i.e. internal or interaction energies
    """
    cols = series.lower().strip().split("/")
    s2a = {"intra":"internal_energies","inter":"interaction_energies"}
    attr = s2a[cols[0]]
    elabel = cols[1]
    etype = cols[2]
  
    label = elabel + "/" + etype + " [kcal/mol]"
    if etype[-1] == "b" :
      eattr = "back"
      etype = etype[:-1]
    elif etype[-1] == "f" :
      eattr = "forw"
      etype = etype[:-1]
    else :
      eattr = "curr"
      
    if hasattr(results,attr) :
      dict = getattr(results,attr)
      if elabel in dict :
        for e in dict[elabel] : 
          if e.type.lower() == etype : return getattr(e,eattr),label
    return None

  def parse_energyresults(series,results) :
    """
    Parse a EnergyResults object
    """
    if series[-1] == "b" :
      eattr = "back"
      etype = series[:-1]
    elif series[-1] == "f" :
      eattr = "forw"
      etype = series[:-1]
    else :
      eattr = "curr"
      etype = series
    label = series + " [kcal/mol]"
    if hasattr(results,etype) :
      return getattr(getattr(results,etype),eattr),label
    else :
      return None
      
  def parse_feenergy(series,results) :
    """
    Parse the feenergy attribute
    """
    if not hasattr(results,"feenergies") : return None
    
    cols = series.lower().strip().split("/")
    if len(cols) == 1 :
      ys  = [results.feenergies[l] for l in sorted(results.feenergies.keys())]
      labels = ["energy at %.3f [kcal/mol]"%l for l in sorted(results.feenergies.keys())]
      return ys,labels
    else :
      ys = [results.feenergies[float(l)] for l in cols[1].split(",") if float(l) in results.feenergies]
      labels = ["energy at %.3f [kcal/mol]"%float(l) for l in cols[1].split(",") if float(l) in results.feenergies]
      if len(ys) == 0 :
        return None
      else :
        return ys,labels  

  # Setup parser for special keywords
  special_parse = {}
  for attr in ["total","capenergy","extraenergy"] :
    special_parse[attr] = parse_energyresults
    special_parse[attr+"f"] = parse_energyresults
    special_parse[attr+"b"] = parse_energyresults
  special_parse["intra"] = parse_compound_series
  special_parse["inter"] = parse_compound_series
  special_parse["feenergy"] = parse_feenergy
  
  # Setup units for special keywords, if not defined here, the unit will be kcal/mol
  special_units = {}
  special_units["solventson"] = ""
  special_units["lambdareplica"] = ""

  ys = []
  labels = []
  for s in args.series :
    resp = None
    attr = s.strip().split("/")[0]
    if attr in special_parse :
      resp = special_parse[attr](s,results)
    else :
      if hasattr(results,s) : 
        if s in special_units :
          l = s + special_units[s]
        else :
          l = s + " [kcal/mol]"
        resp = getattr(results,s),l
    if resp is not None :
      y,label = resp
      if isinstance(y,list) :
        ys.extend(y)
        labels.extend(label)
      else :
        ys.append(y)
        labels.append(label)
    else :
      print "Skipping non-exisiting series %s"%s
  return ys,labels

def _label0(label) :
  """
  Removes the unit from a label
  """
  if label.find("[") == -1 : return label 

  l = "_".join(label.split()[:-1])
  l.replace("/","_")    
  return l

def plot_series(ys,x0s,labels,plotkind,outprefix) :
  """
  Plot a series
  
  Parameters
  ----------
  ys : list of Numpy arrays
    the data series
  x0s : list of int
    the length of the equilibration period for data series
  labels : list of strings
    the labels for the data series
  plotkind : string
    the type of plot to create, can be either sep, sub, single, single_first0 or single_last0
  outprefix : string
    the prefix of the created png-files
  """
  if plotkind not in ["sep","sub","single","single_first0","single_last0"] : plotkind == "sep"
  if len(ys) == 1 or plotkind is None : plotkind = "single" # For a single series, there is only one kind
  ncols = int(np.ceil(len(ys)/2.0))

  ys = np.array(ys)
  x = np.arange(1,ys.shape[1]+1)
  
  # Subtract either first or last point from each series
  if plotkind.startswith("single_") :
    for i in range(ys.shape[0]) :
      if plotkind == "single_first0" :
        ys[i,:] = ys[i,:] - ys[i,0]
      else :
        ys[i,:] = ys[i,:] - ys[i,-1]

  if plotkind != "sep" :
    currfig = plt.figure(1)

  if plotkind != "sub" :
    ymax = ys.max() + 0.1*(ys.max()-ys.min())
    ymin = ys.min() - 0.1*(ys.max()-ys.min())

  for i,(y,x0,label) in enumerate(zip(ys,x0s,labels)) :    
    if plotkind == "sep" :
      currfig = plt.figure(i+1)
    if plotkind == "sub" :
      ax = currfig.add_subplot(2,ncols,i+1)
      ymax = y.max() + 0.1*(y.max()-y.min())
      ymin = y.min() - 0.1*(y.max()-y.min())
    else :
      ax = currfig.gca()

    ax.plot(x,y,label=label,color=simulationobjects.color(i))
    ax.plot([x0,x0],[ymin,ymax],'--',color=simulationobjects.color(i))

    if not plotkind.startswith("single") :
      ax.set_ylim([ymin,ymax])
      ax.set_xlabel("Snapshot")
      ax.set_ylabel(label)
    if plotkind == "sep" :
      currfig.savefig(outprefix+"_"+_label0(label)+".png",format="png")
  
  if plotkind.startswith("single") :
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel("Snapshot")
    if ys.shape[0] == 1 :
      ax.set_ylabel(labels[0])
    else :
      ax.set_ylabel("Data")
      ax.legend()
  if plotkind == "sub" :
    currfig.tight_layout()
  
  if plotkind != "sep" :
    currfig.savefig(outprefix+".png",format="png")
    currfig.show()
    print "\nType enter to quit\n>",
    raw_input()

def write_series(ys,x0s,labels,filekind,outprefix) :
  """
  Write a series to disc
  
  Parameters
  ----------
  ys : list of Numpy arrays
    the data series
  x0s : list of int
    the length of the equilibration period for data series
  labels : list of strings
    the labels for the data series
  filekind : string
    the type of file to write, can be either sep or single
  outprefix : string
    the prefix of the created files
  """
  if len(ys) == 1 or filekind is None : filekind  = "single" # For a single series, there is only one kind
  if filekind not in ["sep","single"] : filekind = "single" # For files the kinds sub, single, single_first0 and single_last0 all means the same

  if filekind == "sep" :
    for y,x0,label in zip(ys,x0s,labels) :      
      with open(outprefix+"_"+_label0(label)+".dat","w") as f :
        f.write("#Data for %s\n"%label)
        f.write("#Equilibration time: %d\n"%x0)
        for i,yi in enumerate(y) :
          f.write("%d %.5f\n"%(i+1,yi))
  else :
    with open(outprefix+".dat","w") as f :
      f.write("#Data for %s\n"%"\t".join(labels))
      f.write("#Equilibration time: %s\n"%"\t".join("%d"%x0 for x0 in x0s))
      for i in range(ys[0].shape[0]) :
        f.write("%d %s\n"%(i+1,"\t".join("%.5f"%y[i] for y in ys)))

#################################
# Wizard routines for user input
#################################

def _select_series(results) :
  """
  Prompts the user for series to plot and write

  Parameters
  ----------
  results : SnapshotResults object
    all the results

  Returns
  -------
  list of string
    the selected series
  """
  selection = []

  print "Select one or several of the following series:"
  print "----------------------------------------------"
  singles = []
  for attr in ["backfe","forwfe","gradient","agradient","lambdareplica","solventson"] :
    if hasattr(results,attr) : singles.append(attr)
  for attr in ["total","capenergy","extraenergy"] :
    if hasattr(results,attr) : singles.append(attr+"[b|f]")
  if len(singles) > 0 : 
    print "Single valued series: "
    for i in range(0,len(singles),5) :
      print "".join("%-15s"%s for s in singles[i:i+5])
    print "(type e.g. total)"
  if hasattr(results,"internal_energies") :
    print "\nInternal energies:"
    for elabel in results.internal_energies :
      print "%-10s : %s"%(elabel,", ".join(e.type.lower() for e in results.internal_energies[elabel]))
    print "(type e.g. intra/protein1/sum)"
  if hasattr(results,"interaction_energies") :
    print "\nInteraction energies:"
    for elabel in results.interaction_energies :
      print "%-20s : %s"%(elabel,", ".join(e.type.lower() for e in results.interaction_energies[elabel]))
    print "(type e.g. inter/solvent-solvent/sum)"
  if hasattr(results,"feenergies") :
    print "\nEnergies at lambda-values: %s"%(", ".join("%.3f"%l for l in sorted(results.feenergies.keys())))        
    print "(type e.g. feenergy or feenergy/0.000,1.000)"
  #if hasattr(series,"thetavals") :
  #  print "\nTheta values for GC-solute: %s"%(", ".join("%d"%i for i in range(1,len(series.thetavals)+1)))        
    
  print "\n> ",
  instr = raw_input()
  while len(instr) > 0 :
    selection.append(instr)
    print "> ",
    instr = raw_input()
  return selection

def _select_plot() :
  """
  Prompt the user to select a plot kind

  Returns 
  -------
  string
    the plot/write kind
  """
  print "\nHow do you want to plot/write the multiple series?"
  print "1) Separate plots (default)"
  print "2) Sub-plots"
  print "3) Single plot"
  print "4) Single plot + subtract first snapshot"
  print "5) Single plot + subtract last snapshot"
  print "\n> ",
  instr = raw_input()
  if instr == "2" : return "sub"
  elif instr == "3" : return "single"
  elif instr == "4" : return "single_first0"
  elif instr == "5" : return "single_last0"
  else : return "sep"

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to analyze and plot a time series")
  parser.add_argument('-f','--file',help="the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-o','--out',help="the prefix of the output figure. Default is results. ",default="results")
  parser.add_argument('-s','--series',nargs="+",help="the series to analyze")
  parser.add_argument('-p','--plot',choices=["sep","sub","single","single_first0","single_last0"],help="the type of plot to generate for several series")
  parser.add_argument('--nperm',type=int,help="if larger than zero, perform a permutation test to determine equilibration, default=0",default=0)
  parser.add_argument('--threshold',type=float,help="the significant level of the equilibration test, default=0.05",default=0.05)
  args = parser.parse_args()
  
  # Read the input file from disc
  results_file = simulationobjects.ResultsFile()
  results_file.read(filename=args.file)
  results = results_file.make_series() # This puts all data into Numpy arrays
  
  # Select which series to plot
  if args.series is None :
    args.series = _select_series(results)
  if len(args.series) == 0 :
    print "No series selected! Exiting"
    quit()
  
  # Parse the user selected series into NumpyArray with the data and an appropriate label
  ys,labels = parse_series(args.series,results)
  if len(ys) == 0 :
    print "No series was parsed from the selection! Exiting"
    quit()
  
  # Compute the equilibration time for each series to plot
  x0s = np.zeros(len(ys))
  x = np.arange(1,ys[0].shape[0]+1)
  print ""
  for i,(y,label) in enumerate(zip(ys,labels)) :
    x0s[i] = find_equilibration(y,x,nperm=args.nperm,threshold=args.threshold)
    print "Equilibration found at snapshot %d for %s"%(x0s[i],_label0(label))

  # Select what kind of plot to make for multiple series
  if len(ys) > 1 and args.plot is None :
    args.plot = _select_plot()
    
  # Plot the series
  plot_series(ys,x0s,labels,args.plot,args.out)

  # Write the series to disc
  write_series(ys,x0s,labels,args.plot,args.out)

