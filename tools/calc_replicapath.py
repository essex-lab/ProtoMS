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

def replica_path(filenames,evalthis,replicakind="lambda") :
  """
  Extract replica path from a set of results files

  Parameters
  ----------
  filenames : list of strings
    the files to read
  evalthis : list
    labels to evaluate the paths on, all files are read
    but paths are only made for the replicas with these labels
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
  else :
    raise simulationobjects.SetupError("Have not implemented replica path analysis for %s"%replicakind)

  paths = []
  labels   = []
  for filename in filenames :
    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=filename)
    if hasattr(results_file.snapshots[0],replica_attr) :
      path = np.zeros(len(results_file.snapshots),int)-1
      label = getattr(results_file.snapshots[0],label_attr)
      if label in evalthis :
        for i,snap in enumerate(results_file.snapshots) :
          path[i] = getattr(snap,replica_attr)
      paths.append(path)
      labels.append(label)

  paths = np.array(paths)
  labeled_paths = np.zeros(paths.shape)
  for i in range(paths.shape[0]) :
    for j in range(paths.shape[1]) :
      if paths[i,j] < 0 : continue
      labeled_paths[i,j] = labels[paths[i,j]-1]
  return labeled_paths,labels

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to analyze and plot a replica paths")
  parser.add_argument('-f','--files',nargs="+",help="the name of the files to analyse")
  parser.add_argument('-p','--plot',type=float,nargs="+",help="the replica values to plot")
  parser.add_argument('-k','--kind',choices=["lambda"],help="the kind of replica to analyze",default="lambda")
  parser.add_argument('-o','--out',help="the prefix of the output figure. Default is replica_path. ",default="replica_path.png")
  args = parser.parse_args()

  # Extract paths and labels from the input files
  paths,labels = replica_path(args.files,args.plot)

  # Plot them and save as a png-file
  x = np.arange(1,paths.shape[1]+1)
  for i,(path,label) in enumerate(zip(paths,labels)) :
    if label in args.plot :
      print label
      plt.plot(x,path,color=simulationobjects.color(i))
  plt.ylim([-0.1,1.1])
  plt.yticks(labels)
  plt.ylabel(args.kind.capitalize())
  plt.xlabel("Snapshot")
  plt.savefig(args.out,format="png")
  
  
