# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to analyze and plot replica paths

This module defines the following public functions:
replica_path

Can be executed from the command line as a stand-alone program
"""

import os
import logging

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pylab as plt

import simulationobjects

logger = logging.getLogger('protoms')

def replica_path(filenames,replicakind="lambda") :
  """
  Extract replica path from a set of results files

  Parameters
  ----------
  filenames : list of strings
    the files to read
  replicakind : string
    the type of replica, at the moment only lambda is allowed

  Returns
  -------
  numpy array
    the labeled replica path
  list
    list of the labels for all input files
  """
  if replicakind == "lambda" :
    replica_attr = "lambdareplica"
    label_attr   = "lam"
  elif replicakind == "temperature" :
    replica_attr = "temperaturereplica"
    label_attr = "temperature"
  elif replicakind == "rest" :
    replica_attr = "temperaturereplica"
    label_attr = "efftemperature"
  elif replicakind == "global" :
    replica_attr = "globalreplica"
    label_attr = None
  elif replicakind == "B" :
    replica_attr = "gcmcreplica"
    label_attr = "bvalue"
  else :
    raise simulationobjects.SetupError("Have not implemented replica path analysis for %s"%replicakind)

  # Build of list of replica ids in each output file
  # also produce a list of labels, e.g. lambda or temperature
  rawpaths = []
  labels   = []
  for ind,filename in enumerate(filenames) :
    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=filename)
    if hasattr(results_file.snapshots[0],replica_attr) :
      path = np.zeros(len(results_file.snapshots),int)-1    
      label = ind+1
      if label_attr : label = getattr(results_file.snapshots[0],label_attr)
      for i,snap in enumerate(results_file.snapshots) :
        path[i] = getattr(snap,replica_attr)
      rawpaths.append(path)
      labels.append(label)

  # Now find the path for each lambda/temperature value
  rawpaths = np.array(rawpaths,dtype=int)
  labeled_paths = [[] for l in labels]
  # Loop over each possible lambda/temperature value, i.e. label id
  for i in range(rawpaths.shape[0]) :
    # We will try to fill up an array that is as long as the simulation
    while len(labeled_paths[i]) < rawpaths.shape[1] :
      k = len(labeled_paths[i]) # This is the current point in time we will be looking for replica with id = i 
      # Look in the list of each output file as built above
      for j,path in enumerate(rawpaths) :     
        # If the replica id at the current time, k is equal to the label id we are on (i),
        # add the label of file j to the list and break
        if path[k] == i + 1 : 
          labeled_paths[i].append(labels[j])
          break
  return np.asarray(labeled_paths),labels

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to analyze and plot a replica paths")
  parser.add_argument('-f','--files',nargs="+",help="the name of the files to analyse")
  parser.add_argument('-p','--plot',type=float,nargs="+",help="the replica values to plot")
  parser.add_argument('-k','--kind',choices=["lambda","temperature","rest","global", "B"],help="the kind of replica to analyze",default="lambda")
  parser.add_argument('-o','--out',help="the prefix of the output figure. Default is replica_path. ",default="replica_path.png")
  return parser

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  args = get_arg_parser().parse_args()

  # Extract paths and labels from the input files
  paths,labels = replica_path(args.files,args.kind)

  # Plot them and save as a png-file
  x = np.arange(1,paths.shape[1]+1)
  for i,(path,label) in enumerate(zip(paths,labels)) :
    if label in args.plot :
      plt.plot(x,path,color=simulationobjects.color(i))
  psorted = sorted(args.plot)
  prange = max(psorted) - min(psorted)
  plt.ylim([min(psorted)-0.1*prange,max(psorted)+0.1*prange])
  plt.yticks(labels)
  plt.ylabel(args.kind.capitalize())
  plt.xlabel("Snapshot")
  plt.savefig(args.out,format="png")
  
  
