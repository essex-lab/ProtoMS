# Author: Gregory Ross

""" 
Program to analyze GCMC simulations
"""

import logging
import glob
import sys
import os
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit
import simulationobjects
from simulationobjects import ResultsFile
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pyplot as plt

logger = logging.getLogger('protoms')

def logistic_single(B,dF):
    kT = 0.5924837
    return 1/(1+np.exp(dF/kT - B)) 

def fit_single(B,N,initial):
     def errfunc(x): return np.sum((N  - logistic_single(B,x) )**2)
     solution = optimize.minimize(errfunc, x0 =initial,method="BFGS") 		# BFGS is faster than "Powell"
     return solution.x

# A Gaussian smoother.
def gaussian_smooth(x,y,sigma=None):
    if sigma is None:
        sigma = 0.05
    smooth_y = np.zeros(len(y))
    for i in range(len(x)):
        weights = np.exp(-((x[i] - x)**2)/(2*sigma))
        weights = weights/np.sum(weights)
        smooth_y[i] = np.sum(weights*y)
    return smooth_y

def plot_FitPercentiles(B,N,dF_samples,resolution=None,level1=None,level1_colour=None,level2=None,level2_colour=None,median_colour=None,smoothness=None):
    # Setting the defaults:
    if resolution == None: resolution = 50					# Number of points with which the lines will be drawn.
    if level1 == None: level1=0.50							# Default is also shade the region that contains 50% of the data
    if level2 == None: level2=0.90				        	# Default is NOT to shade a region at a lower percentile level.
    if level1_colour == None: level1_colour="orange"		# Default is also shade the region that contains 50% of the data
    if level2_colour == None: level2_colour="gray"			# Default is to also shade the region that contains 95% of the data
    if median_colour == None: median_colour="red"			# Colour of the median
    if smoothness == None: smoothness = 0.05				# Degree to which the lines will be smoothed (for aesthetic purposes only).
    # Generating the predictions of the different gcmc_models:
    yvals = np.zeros((resolution,len(dF_samples)))						# Will hold all the predictions from the different models
    x = np.linspace(start=B.min(),stop=B.max(),num=resolution)		# The x values that will be plotted.
    for i in range(len(dF_samples)):
        yvals[:,i] = logistic_single(x,dF_samples[i])
    # Generating and smoothing the confidense intervals.
    low, high = (1-level1)/2, (1-level1)/2 + level1
    y_median, y_low, y_high = np.percentile(yvals,50,axis=1),np.percentile(yvals,100*low,axis=1),np.percentile(yvals,100*high,axis=1)
    low, high = (1-level2)/2, (1-level2)/2 + level2
    y_median, y_lowest, y_highest = np.percentile(yvals,50,axis=1),np.percentile(yvals,100*low,axis=1),np.percentile(yvals,100*high,axis=1)
    smooth_median =gaussian_smooth(x,y_median,sigma=smoothness*3)		# The median is made smoother than the error bars.
    smooth_low = gaussian_smooth(x,y_low,sigma=smoothness)
    smooth_high = gaussian_smooth(x,y_high,sigma=smoothness)
    smooth_lowest = gaussian_smooth(x,y_lowest,sigma=smoothness)
    smooth_highest = gaussian_smooth(x,y_highest,sigma=smoothness)
    # Plotting:
    space = 0.07			# The amount of white space to the side of the x and y axis as a fraction of the ranges.
    plt.fill_between(x,smooth_lowest,smooth_highest, facecolor=level2_colour,linewidth=0,alpha=0.3,interpolate=True)
    plt.fill_between(x,smooth_low,smooth_high, facecolor="white",linewidth=0,interpolate=True)
    plt.fill_between(x,smooth_low,smooth_high, facecolor=level1_colour,linewidth=0,alpha=0.4,interpolate=True)
    plt.xlim(B.min()-space*np.ptp(B),B.max()+space*np.ptp(B))
    plt.ylim(N.min()-space*np.ptp(N),N.max()+space*np.ptp(N))
    plt.scatter(B, N,color="black",s=40)
    plt.plot(x, smooth_median,color=median_colour,linewidth=3)
    return plt

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate coupling free energy GCMC simulations of a single water molecule using the logistic fitting method by G. A. Ross et al.")
  parser.add_argument('-d','--directories',nargs="+",help="the directories containing the GCMC simulation data, default=None",default=None)
  parser.add_argument('-f','--file',help="the name of the file to analyse, default=results",default="results")
  parser.add_argument('-s','--skip',help="the number of initial snapshots that will be discarded, default=0",type=int,default=0)
  parser.add_argument('-c','--calc',help="whether to calculate the coupling free energy by fitting a logistic curve to titration data, default=True",type=bool,default=True)
  parser.add_argument('-p','--plot',action='store_true',help="whether to plot the titration data and fitted curve, default is no plotting",default=False)
  parser.add_argument('--guess',help="initial guess of the coupling free energy. Use to refine logisitic fitting, default=-6.2 kcal/mol",type=float,default=-6.2)
  parser.add_argument('--excess',action='store_true',help="calculate the average excess chemical potential between 2 Adams values, default=False",default=False)
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
        resultsfiles = glob.glob(folder+ "/"+args.file+"*")
        if len(resultsfiles) > 1: # It is assumed the results are from ProtoMS 2.
            results.read([folder,args.file])
        elif len(resultsfiles)==1: # It is assumed the results are from ProtoMS 3.
            results.read(folder+ "/"+args.file)
        else:
           print "\nError. No results file matching %s. Exiting program\n" % folder+"/"+args.file 
           sys.exit()
        solventson = np.array([snap.solventson for snap in results.snapshots])# Counting how many molecules have been inserted.
        mean_on = np.mean(solventson[args.skip:solventson.size])
        N.append(mean_on)
        adams = np.array([snap.bvalue for snap in results.snapshots])	# Recording the B-value for each window (should be constant, but averaging just in case).
        mean_adams = np.mean(adams[(args.skip-1):adams.size])
        B.append(mean_adams)
  print "\n...GCMC data has been read."

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

  # Fitting to the data:
  if args.calc==True:
    dF = fit_single(B_sorted,N_sorted,args.guess)
    # Bootstrap fitting to obtain error estimate:
    num_boots = 1000
    dF_boots = np.zeros(num_boots)    
    indices = range(B.size)
    for boot in range(num_boots): 
      sample_inds = np.random.choice(indices,size=B.size)
      B_sample = B_sorted[sample_inds]
      N_sample = N_sorted[sample_inds]
      dF_boots[boot] = fit_single(B_sample,N_sample,dF)
    print "\nFREE ENERGY ESTIMATES:"
    print "  Least squares estimate = %.2f kcal/mol" % dF
    print "  Bootstrap estimate     = %.2f +/- %.2f kcal/mol\n" % (np.mean(dF_boots),np.std(dF_boots))  
    print "The free energy stated is the free energy to transfer from ideal gas to the GCMC region. The binding free energy is the difference between this value and the hydration free energy of the molecule."

  # Plotting
  if args.plot==True:
    plt.figure("Least squares fit")
    plt.scatter(B_sorted, N_sorted,color="k",s=40)
    if args.calc==True: 
      Bx = np.linspace(start=B.min(),stop=B.max(),num=500)
      plt.plot(Bx,logistic_single(Bx,dF),color="red",linewidth=3)
    plt.xlabel("Adams parameter (B)",fontsize=15)
    plt.ylabel("Average amount of water",fontsize=15)
    plt.suptitle("Fitted free energy and logistic curve",fontweight="bold")
    plt.show(block=False)

    plt.figure("Bootstrap Samples")
    currfig = plot_FitPercentiles(B,N,dF_boots,level1=0.50,level2=0.90,smoothness=0.01)
    currfig.xlabel("Adams parameter (B)",fontsize=15)
    currfig.ylabel("Average amount of water",fontsize=15)
    currfig.suptitle("Median and percentiles of fitted free energy",fontweight="bold")
    currfig.show(block=False)

  # Using consistent values of dG to average:
  if args.excess==True: 
    notzero = np.where(N_sorted>1E-6)
    excess = (B_sorted[notzero] - np.log(N_sorted[notzero]))*0.592
    plt.figure("Excess chemical potential for each Adams value")
    plt.plot(B_sorted[notzero], excess,color="k",linewidth=3)
    plt.xlabel("Adams parameter (B)",fontsize=15)
    plt.ylabel("Excess chemical potential (kcal/mol)",fontsize=15)
    plt.show(block=False)
    print "\nPlease enter the lower and upper Adams values between which the mean excess chemical potential will be calculated. Where the excess chemical potential appears flat, it can be used to approximate the free energy to transfer the molecule from ideal gas to the GCMC region.\n>"
    inputrange = raw_input()
    inputrange = inputrange.split(" ")
    if float(inputrange[0]) < float(inputrange[1]):
      B_low = float(inputrange[0])
      B_high = float(inputrange[1])
    else:
      B_low = float(inputrange[1])
      B_high = float(inputrange[0])
    indices = (B_sorted[notzero] >= float(B_low))*(B_sorted[notzero] <= float(B_high)) # The indices from which the mean excess chemical potential will be calculated.
    print "\nAverage excess chemical potential of region = %.2f kcal/mol, with standard error = %.2f kcal/mol. Number of data points used = %i.\n" % (np.mean(excess[indices]), np.std(excess[indices])/np.sqrt(np.sum(indices)), np.sum(indices) )

  if args.plot==True or args.excess==True:
    print "\nType enter to quit\n>"
    raw_input()

