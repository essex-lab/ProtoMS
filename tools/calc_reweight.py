# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to reweight replica exchange results
to use information from all replicas.

It can be excuted from the command line as
a stand alone program.
"""

import numpy as np
import simulationobjects
import math

import os
import logging

logger = logging.getLogger('protoms')


def replica_histogram(energy="total",replica="temperature",variable="theta",ewidth=100.0,rwidth=1.0,results=["results"],skip=0) :
  """
  Returns a numpy matrix which includes effectively the data
  of a 2D histogram in the "energy" and "replica" dimensions,
  as well as the average value of the "variable" for each bin.

  Parameters
  ----------
  energy: string
    the kind of energy that will be stored
  replica: string
    the kind of replica that will be plotted
  variable: string
    the variable which average value per bin
    will be stored
  ewidth: float
    width of the energy bins
  results: string
    path of the results files
  skip: integer
    number of snapshots to skip

  Returns
  -------
  2D numpy histogram
    the histogram over energy and replica
  numpy 1d array
    the mean of the variable for each energy bin
  """

  results_files = []

  for eachres in results :
    each_file = simulationobjects.ResultsFile()
    each_file.read(filename=eachres,skip=skip)
    results_files.append(each_file)

#   This might be the most general way, but I am not sure if this is the best one, because I could directly store thetas in terms of temperatures. Same for temperatures.

  arrenergies = []
  arrvariable = []
  arrreplicas = []

  for results_file in results_files :
    for snap in results_file.snapshots :
      if energy == "total" :
        arrenergies.append(snap.total.curr)
      else :
        simulatioobjects.SetupError("This tipe of energy is not supported.")

#   I also need to know which theta I want among all those that are there! (for the different solutes)!!!!!!

      if variable == "theta" :
        if snap.thetasolvals :
          arrvariable.append(snap.thetasolvals[0])
        else :
          arrvariable.append(snap.thetavals[0])
      else :
        simulatioobjects.SetupError("This tipe of variable is not supported.")

      if replica == "temperature" :
        arrreplicas.append(snap.temperature)
      else :
        simulatioobjects.SetupError("This tipe of replica is not supported.")

  arrenergies = np.array(arrenergies)
  arrvariable = np.array(arrvariable)
  arrreplicas = np.array(arrreplicas)

  print "arrenergies"
  print arrenergies

  nebins = int(math.ceil((np.amax(arrenergies)-np.amin(arrenergies)+1)/ewidth))
  nrbins = int(math.ceil((np.amax(arrreplicas)-np.amin(arrreplicas)+1)/rwidth))

  print nebins
  print nrbins

  hist_enre = np.histogram2d(arrenergies,arrreplicas,bins=[nebins,nrbins])

  varmeans = np.zeros((nebins))

  for i in range(nebins) :

    varmeans[i] = np.mean(np.array([float(var) for j,var in enumerate(arrvariable) if hist_enre[1][i] <= arrenergies[j] < hist_enre[1][i+1] ]))

  return hist_enre,varmeans




def reweighted_variable(varmeans,hist,replica="temperature",feguess=10.0,tolerance=0.5,gcorr=1.0,iterations=500,main=None) :
  """
  A function to calculate the ensemble average of a variable
  following the WHAM reweightening from a 2d histogram of 
  energies and some bias or replicas and the mean of the
  variable for each energy bin.

  Parameters
  ----------
  varmeans: 1d numpy array
    the mean of the variable for each energy bin
  hist: 2d numpy histogram
    the histogram over energy and replicas / bias
  replica: string
    the type of replica over which WHAM will be applied
  feguess: float
    initial guess of the free energy associated
    with each replica / bias bin
  tolerance: float
    the maximum difference between guessed free energy
    and that calculated from the density of states
  gcorr: float
    correction that accounts for statistical inefficiency
  iterations: integer
    maximum number of iterations trying to reach tolerance
  main: string or None
    the replica towards which the reweightening is applied

  Returns
  -------
  float
    reweighted ensemble average of the variable value
  """

  histval = hist[0]

  fe = np.array([np.zeros(histval.shape[1])+feguess])

# The line below is specific for temperatures!!

  kb = 0.0019872041
  C2K = 273

  ewidth = (hist[1][1]-hist[1][0])
  rwidth = (hist[2][1]-hist[2][0])

  if replica == "temperature" :

    beta = np.array([1/(kb*(np.delete(hist[2],-1)+rwidth/2+C2K))])

  else :
    
    beta = np.delete(hist[2],-1)+rwidth/2


  energy = np.swapaxes(np.array([np.delete(hist[1],-1)+ewidth/2]),0,1)

# A constant (C) is needed to apply both in the calculation of free energies and density of states, to avoid the exponential of a huge negative number

  C = (hist[1][0]+hist[1][-1])/2

  print "beta"
  print beta

  print "ENERGY HIST"
  print hist[1]
  print hist[0]

  iteration = 0

  fe_new = fe

  for iteration in range(iterations) :

    fe = fe_new

    density_est = np.swapaxes(np.array([np.sum(gcorr*histval,axis=1)/np.sum(gcorr*np.array([np.sum(histval,axis=0)])*ewidth*np.exp(fe-beta*(energy-C)),axis=1)]),0,1)

    print 

    fe_new = (-1)*np.array([np.log(np.sum(density_est*ewidth*np.exp((-1)*beta*(energy-C)),axis=0))])

    print energy-C

    print "CHECK BELOW"

    print np.exp(fe-beta*(energy-C))

    if iteration%(iterations/10) == 0 :
      logger.info("------------")
      logger.info("Iteration %d:"%(iteration+1))
      logger.info("Free energies: ")
      logger.info(fe_new.reshape(-1))
      logger.info("Desnities of states: ")
      logger.info(density_est.reshape(-1))

    print "BAS BELOW"

    print np.absolute(fe-fe_new)

    if (np.absolute(fe-fe_new)<tolerance).all() :
      break

  print "density_est"
  print density_est

  print "free energy"
  print fe_new

  if isinstance(main,str) :

    try :
      mybeta = beta[0][int(main)-1]
      logger.info("Replica %s chosen as main."%main)
    except:
      mybeta = beta[0][0]
      logger.info("Replica %s could not be interpreted. Lowest replica chosen as main."%main)

  else :
    mybeta = beta[0][0]
    

  vararr = np.swapaxes(np.array([varmeans]),0,1)

  enav_var = np.sum(density_est*ewidth*np.exp(mybeta*energy)*vararr,axis=0)/np.sum(density_est*ewidth*np.exp(mybeta*energy),axis=0)

  print "Variable reweighted:"
  print enav_var

  return enav_var
  



#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  parser = argparse.ArgumentParser(description="Program to reweight energies from different replicas with WHAM.")
  parser.add_argument('-f','--files',nargs="+",help="the path of the files to analyse. Default is results. ",default=["results"])
  parser.add_argument('-r','--replicas',choices=["temperature"],help="which replicas we are using",default="temperature")
  parser.add_argument('-v','--variable',choices=["theta"],help="the variable to reweight",default="theta")
  parser.add_argument('-s','--skip',help="number of snapshots to skip from the results file. Default=0.",default=0)
  parser.add_argument('-eb','--energybin',help="the width of the bins for the energy",default=100.0)
  parser.add_argument('-rb','--replicabin',help="the width of the bins for the temperature",default=1.0)
  parser.add_argument('--iterations',help="the number of WHAM iterations allowed. Default=500",default=500)
  parser.add_argument('--mainreplica',help="the replica of interest. By default, the lowest one is chosen.",default=None)
  args = parser.parse_args()

  try :
    skip = int(args.skip)
  except :
    simulationobjects.SetupError("Expected integer after flag '-s' / '--skip'.")

  hist_enre, varmeans = replica_histogram(replica=args.replicas,variable=args.variable,ewidth=args.energybin,rwidth=float(args.replicabin),results=args.files,skip=skip)

  reweighted_variable(varmeans, hist_enre,iterations=args.iterations,main=args.mainreplica)












        
