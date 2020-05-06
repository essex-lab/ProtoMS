# Authors: -The ProtoMS developers team-
#           Ana I. Cabedo Martinez
#           Gregory A. Ross


"""
Program to reweight replica exchange results
to use information from all replicas.

It can be excuted from the command line as
a stand alone program.
"""

import matplotlib.pyplot as plt
import numpy as np
import simulationobjects
import math
import sys

import os
import logging

logger = logging.getLogger('protoms')



def replica_histogram(energy="total",replica="temperature",variable="theta",which='1',ewidth=100.0,rwidth=1.0,results=["results"],skip=0,plothist=True) :
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
  which: string
    which of the instances of variables is to
    be stored, if more than one is possible
  ewidth: float
    width of the energy bins
  results: string
    path of the results files
  skip: integer
    number of snapshots to skip
  plothist: boolean
    whether to plot the 2dhistogram

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
        raise simulatioobjects.SetupError("This tipe of energy is not supported.")

#   I also need to know which theta I want among all those that are there! (for the different solutes)!!!!!!

      if variable == "theta" :
        if snap.thetasolvals :
          arrvariable.append(snap.thetasolvals[int(which)-1])
        else :
          arrvariable.append(snap.thetavals[int(which)-1])
      elif variable == "totalenergy" :
          arrvariable = arrenergies
      else :
        raise simulatioobjects.SetupError("This tipe of variable is not supported.")

      if replica == "temperature" :
        arrreplicas.append(snap.temperature)
      else :
        raise simulatioobjects.SetupError("This tipe of replica is not supported.")

  arrenergies = np.array(arrenergies)
  arrvariable = np.array(arrvariable)
  arrreplicas = np.array(arrreplicas)

  logger.debug("Maximum energy found:")
  logger.debug(arrenergies.max())

  logger.debug("Minimum energy found:")
  logger.debug(arrenergies.min())

  replicas = np.unique(arrreplicas)

  nebins = int(math.ceil((np.amax(arrenergies)-np.amin(arrenergies))/ewidth))

  # if no width for replicas is provided, the minimum distance between any of them is taken
  if not rwidth :
    try : 
      rwidth = math.fabs(np.subtract(replicas[:-1],replicas[1:]).min())
    except :
      raise simulationobjects.SetupError("Automatic replica bin width could not be stated. Please, provide an argument for --replicabin.")
    
  arr_rbin_edges = np.array([replica-(rwidth/2) for replica in replicas] + [replicas[-1]+(rwidth/2)])
  arr_ebin_edges = np.array([np.amin(arrenergies)+ewidth*i for i in range(nebins+1)])

  hist_enre = np.histogram2d(arrenergies,arrreplicas,bins=[arr_ebin_edges,arr_rbin_edges])

  varmeans = np.zeros((nebins))

  for i in range(nebins) :

    varmeans[i] = np.mean(np.array([float(var) for j,var in enumerate(arrvariable) 
                  if (hist_enre[1][i] <= arrenergies[j] < hist_enre[1][i+1] and i != nebins-1) or (hist_enre[1][i] <= arrenergies[j] <= hist_enre[1][i+1] and i == nebins-1) ]))


  if plothist :
    filename = "histogram.png"
    fig = plt.figure()
    H = np.rot90(hist_enre[0])
    H = np.flipud(H)
    plt.pcolormesh(hist_enre[1],hist_enre[2],H)
    plt.xlabel('energy (kcal/mol)')
    plt.ylabel('temperature (Celsius)')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(filename)
    logger.info("Histogram saved in '%s'"%filename)

  return hist_enre,varmeans



def set_variables(hist,replica,main) :
  """
  A function ot set the variables required
  by both twham and weight functions

  Parameters
  ----------
  hist: 2d numpy histogram
    the histogram over energy and replicas / bias
  replica: string
    the type of replica over which WHAM will be applied
  main: string or None
    the replica towards which the reweightening is applied

  Returns
  -------
  array
    values of the histogram in hist
  array
    middle point of each energy bin
  array
    value of replica bins, transformed if required
  float
    value of the replica of interest, transformed if required
  float
    width of energy bins
  float
    width of replica bins
  """
  
  histval = hist[0]

  ewidth = (hist[1][1]-hist[1][0])
  rwidth = (hist[2][1]-hist[2][0])

  energy = np.expand_dims(np.delete(hist[1],-1)+ewidth/2,1)

  if replica == "temperature" :

    kb = 0.0019872041
    C2K = 273 

    beta = np.expand_dims(1/(kb*(np.delete(hist[2],-1)+rwidth/2+C2K)),0)

  else :
    
    beta = np.expand_dims(np.delete(hist[2],-1)+rwidth/2,0)

  if isinstance(main,str) :

    try :
      mybeta = float(main)
      if replica == "temperature" : mybeta = 1/(kb*(float(main)+C2K))
      logger.info("Replica %s chosen as main."%main)
    except:
      mybeta = beta[0][0]
      logger.info("Replica %s could not be interpreted. Lowest replica chosen as main."%main)

  else :
    mybeta = beta[0][0]

  return histval,energy,beta,mybeta,ewidth,rwidth



def twham (hist,replica="temperature",tolerance=0.001,gcorr=1.0,iterations=500,out_freq=1000,main=None) :
  """
  A function to calculate the WHAM reweightening from a
  2d histogram of energies and some bias or replicas.

  Parameters
  ----------
  hist: 2d numpy histogram
    the histogram over energy and replicas / bias
  replica: string
    the type of replica over which WHAM will be applied
  tolerance: float
    the maximum difference between guessed free energy
    and that calculated from the density of states
  gcorr: float
    correction that accounts for statistical inefficiency
  iterations: integer
    maximum number of iterations trying to reach tolerance
  out_freq: integer, optional
    the frequency of writting output from the iterations
  main: string or None
    the replica towards which the reweightening is applied

  Returns
  -------
  array
    weights
  array
    middle point of each energy bin
  array
    value of replica bins, transformed if required
  float
    value of the replica of interest, transformed if required
  """

  histval,energy,beta,mybeta,ewidth,rwidth = set_variables(hist,replica,main)

  # Estimated free energy guess (feguess) to avoid numerical inestabilities when calculation density of states and free energy

  feguess =  beta*(hist[1][0]+hist[1][-1])/2

  fe = feguess

  iteration = 0

  fe_new = fe

  Nk = np.expand_dims(np.sum(histval,axis=0),0)

  for iteration in range(iterations) :

    fe = fe_new

    density_est = np.expand_dims(np.sum(histval/gcorr,axis=1)/np.sum(Nk*ewidth*np.exp(fe-beta*energy)/gcorr,axis=1),1)

    fe_new = np.expand_dims((-1)*np.log(np.sum(density_est*ewidth*np.exp((-1)*beta*energy),axis=0)),0)

    if iterations < out_freq or iteration%(iterations/out_freq) == 0 :
      logger.debug("------------")
      logger.debug("Iteration %d:"%(iteration+1))
      logger.debug("Free energies: ")
      logger.debug(fe_new.reshape(-1))
      logger.debug("Densities of states: ")
      logger.debug(density_est.reshape(-1))

    if (np.absolute(fe-fe_new)<tolerance).all() :
      density_est = np.expand_dims(np.sum(histval/gcorr,axis=1)/np.sum(Nk*ewidth*np.exp(fe_new-beta*energy)/gcorr,axis=1),1)
      break

  # A constant (C) is applied to avoid inestabilities in the exponential. It will cancel in 'weighted_variable'

  C = mybeta*(hist[1][0]+hist[1][-1])/2

  weight = density_est*np.exp((-1)*mybeta*energy+C)

  return weight,energy,beta,mybeta



def weight(hist,replica="temperature",main=None) :
  """
  A function to calculate the WHAM reweightening from a
  2d histogram of energies and some bias or replicas.

  Parameters
  ----------
  hist: 2d numpy histogram
    the histogram over energy and replicas / bias
  replica: string
    the type of replica over which WHAM will be applied
  main: string or None
    the replica towards which the reweightening is applied

  Returns
  -------
  array
    weights
  array
    middle point of each energy bin
  array
    value of replica bins, transformed if required
  float
    value of the replica of interest, transformed if required
  """

  histval,energy,beta,mybeta,ewidth,rwidth = set_variables(hist,replica,main)

  if main is None : main = hist[2][0] - rwidth/2 

  if not isinstance(main,float) : main = float(main)

  myreplica = [ ind for ind,leftedge in enumerate(np.delete(hist[2],-1)) if leftedge < main < hist[2][ind+1] ][0]

  myhistvals = np.expand_dims(histval[:,myreplica],1)

  return myhistvals,energy,beta,mybeta


def weighted_variable(varmeans,weight,main,energy,beta,mybeta,replica) :
  """
  A function to reweight average variables.

  Parameters
  ----------
  varmeans: 1d numpy array
    average variable for each bin
  weight: numpy array
    weight to be applied to each bin in varmeans
  main: string
    the replica towards which the reweightening is applied
  energy: array
    middle point of each energy bin
  beta: array
    value of replica bins, transformed if required
  mybeta: float
    value of the replica of interest, transformed if required
  replica: string
    the type of replica over which WHAM will be applied

  Returns
  -------
  float
    reweighted ensemble average of the variable value
  """


# Curating the values of my variable. If my density of states is 0, the "var" value for this bin will not have an influence in the final weighted average anyway.
# The lines below are hence should not have a mathematical effect and sort computational problems

  for ind,eachweight in enumerate(weight.reshape(-1)) :
    if eachweight == 0 :
      varmeans[ind] = 1.0

  vararr = np.expand_dims(varmeans,1)

  enav_var = np.sum(weight*vararr,axis=0)/np.sum(weight,axis=0)

  logger.info("\nWeighted variable: %.5f"%enav_var)

  return enav_var


#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  parser = argparse.ArgumentParser(description="Program to reweight energies from different replicas with WHAM.")
  parser.add_argument('-f','--files',nargs="+",help="the path of the files to analyse. Default is results. ",default=["results"])
  parser.add_argument('-r','--replicas',choices=["temperature"],help="which replicas we are using",default="temperature")
  parser.add_argument('-v','--variable',choices=["theta","totalenergy"],help="the variable to reweight",default="theta")
  parser.add_argument('-w','--weight',choices=["single","twham"],help="estimate of the density of states (weighting factor). 'single' uses mainreplica, 'twham' uses all replicas. Default=single",default="single")
  parser.add_argument('-s','--skip',help="number of snapshots to skip from the results file. Default=0.",default='0')
  parser.add_argument('-eb','--energybin',help="the width of the bins for the energy",default=100.0)
  parser.add_argument('-rb','--replicabin',help="the width of the bins for the temperature",default=None)
  parser.add_argument('-wv','--whichvar',help="if more than one instance of the variable is in files, which of them is to be reweighted. Default=1",default=1)
  parser.add_argument('--iterations',help="the number of WHAM iterations allowed. Default=500",default=500)
  parser.add_argument('--mainreplica',help="the replica of interest. By default, the lowest one is chosen.",default=None)
  parser.add_argument('--plothistogram',action="store_true",help="whether to plot a histogram of the variable counts per energy and replica bin. Default=false",default=False)
  args = parser.parse_args()

  logger = simulationobjects.setup_logger("calc_reweight_py.log")
  logger.info("")

  try :
    skip = int(args.skip)
  except :
    raise simulationobjects.SetupError("Expected integer after flag '-s' / '--skip'.")

  try :
    iterations = int(args.iterations)
  except :
    raise simulationobjects.SetupError("Value %s for --iterations could not be interpreted as integer"%args.iterations) 

  if args.replicabin != None  :
    try :
      args.replicabin = float(args.replicabin)
    except :
      raise simulationobjects.SetupError("Value %s for --replicabin could not be interpreted as float"%args.replicabin)

  hist_enre, varmeans = replica_histogram(replica=args.replicas,variable=args.variable,which=args.whichvar,ewidth=float(args.energybin),rwidth=args.replicabin,results=args.files,skip=skip)

  if args.weight == "twham" :
    weight,energy,beta,mybeta = twham(hist_enre,replica=args.replicas,iterations=iterations,main=args.mainreplica)
  elif args.weight == "single" :
    weight,energy,beta,mybeta = weight(hist_enre,replica=args.replicas,main=args.mainreplica)
  else :
    raise simulationobjects.SetupError("Value %s for --weight could not be interpreted."%args.weight)

  result = weighted_variable(varmeans,weight,args.mainreplica,energy,beta,mybeta,args.replicas)

  logger.info("")












        
