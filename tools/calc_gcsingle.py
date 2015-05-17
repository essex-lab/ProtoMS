# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Program to analyze GCMC simulations
"""

import logging
import glob
import sys
import os

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import simulationobjects
from simulationobjects import ResultsFile

logger = logging.getLogger('protoms')

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to analyze and plot free energies from GCMC simulations")
  parser.add_argument('-d','--directories',nargs="+",help="the directories containing the GCMC simulation data, default=None",default=None)
  parser.add_argument('-f','--file',help="the name of the file to analyse. Default is results.",default="results")
  parser.add_argument('-s','--skip',help="the number of initial snapshots that will be discarded",type=int,default=0)
  parser.add_argument('-p','--plot',action='store_false',help="whether to plot the estimated excess chemical potential against Adams value",default=True)
  parser.add_argument('-r','--range',nargs=2,help="range of Adams value from which the free energy will be calculated.",default=None)
  args = parser.parse_args()

  dG_hyd = -6.2

  # Read in the GCMC results data from multiple GCMC output folders and calculating the mean number of "on" waters, after skipping a specified number of frames.
  directories = args.directories
  N = []
  B = []
  print "\nREADING GCMC DATA:"
  for dirs in directories:
    folders =  glob.glob(dirs)
    if len(folders)==0:
        print "\nError. No folder(s) matching '%s'. Exiting program.\n" % args.directories
        sys.exit()
    for folder in folders:
        results = ResultsFile()
        resultsfiles = glob.glob(os.path.join(folder,args.file+"*"))
        if len(resultsfiles) > 1: # It is assumed the results are from ProtoMS 2.
            results.read([folder,args.file])
        elif len(resultsfiles)==1: # It is assumed the results are from ProtoMS 3.
            results.read(os.path.join(folder,args.file))
        else:
           print "\nError. No results file matching %s. Exiting program\n" % folder+os.path.sep+args.file 
           sys.exit()
        solventson = np.array([snap.solventson for snap in results.snapshots])# Counting how many molecules have been inserted.
        mean_on = np.mean(solventson[args.skip:solventson.size])
        N.append(mean_on)
        adams = np.array([snap.bvalue for snap in results.snapshots])	# Recording the B-value for each window (should be constant, but averaging just in case).
        mean_adams = np.mean(adams[(args.skip-1):adams.size])
        B.append(mean_adams)

  # Checking to make sure the data makes sense.
  N = np.array(N)
  B = np.array(B)
  if B.size != N.size:
      print "\nFatal error: Number of Adams values must equal mean number of inserted waters.\n"
      sys.exit()

  # Sorting the data
  order = np.argsort(B)
  B_sorted = B[order]
  N_sorted = N[order]
  notzero = np.where(N_sorted>1E-6)
  dG = (B_sorted[notzero] - np.log(N_sorted[notzero]))*0.592

  # Plotting
  if args.plot is not None:
    plt.plot(B_sorted[notzero], dG,color="k",linewidth=3)
    plt.xlabel("Adam's parameter (B)",fontsize=15)
    plt.ylabel("Excess chemical potential (kcal/mol)",fontsize=15)
    plt.show(block=False)

  # Using consistent values of dG to average:
  if args.range is not None:
    if len(args.range) == 2:
      if float(args.range[0]) < float(args.range[1]):
        B_low = float(args.range[0])
        B_high = float(args.range[1])
      else:
        B_low = float(args.range[1])
        B_high = float(args.range[0])
      indices = (B_sorted[notzero] >= float(B_low))*(B_sorted[notzero] <= float(B_high)) # The indices from which the mean excess chemical potential will be calculated.
    else:
      print "Two arguements expected for 'range'. Please re-enter the lower and upper B values between which the average will be calculated."
      quit()
    indices = (B_sorted[notzero] >= float(B_low))*(B_sorted[notzero] <= float(B_high))
  else:
    print "\nEnter the lower and upper Adams values between which the mean excess chemical potential will be calculated. \n>"
    inputrange = raw_input()
    inputrange = inputrange.split(" ")
    if float(inputrange[0]) < float(inputrange[1]):
      B_low = float(inputrange[0])
      B_high = float(inputrange[1])
    else:
      B_low = float(inputrange[1])
      B_high = float(inputrange[0])
    indices = (B_sorted[notzero] >= float(B_low))*(B_sorted[notzero] <= float(B_high)) # The indices from which the mean excess chemical potential will be calculated.


  print "\nEstimated free energy to transfer water from ideal gas = %.2f kcal/mol, with standard deviation = %.2f kcal/mol. Number of data points used = %i.\n" % (np.mean(dG[indices]), np.std(dG[indices]), np.sum(indices) )

  print "\nType enter to quit\n>"
  raw_input()

