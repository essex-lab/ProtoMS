# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to do thermodynamic integration

This module defines a single public function:
ti

Can be executed from the command line as a stand-alone program
"""

import glob
import os

import matplotlib.pyplot as pl
import numpy as np

import simulationobjects

def _print_ene(lam,ene,std,print_uncert,print_lam) :
  """ 
  Print energy row to standard output
  
  Parameters
  ----------
  lam : float
    the lambda value
  ene : float
    the energy value
  std : float
    the uncertainty
  print_uncert : boolean
    if to print the uncertainty
  print_lam : boolean
    if to print the lambda value
  """
  if print_lam : print "%-6.2f"%lam,
  print " %12.4f"%ene,
  if print_uncert : print " %8.4f"%std,
  print ""

def _print_head(lam,ene,std,print_uncert,print_lam) :
  """ 
  Print energy header to standard output

  Parameters
  ----------
  lam : string
    the header of the lambda value
  ene : string
    the header of the energy value
  std : string
    the header of the uncertainty
  print_uncert : boolean
    if to print the uncertainty
  print_lam : boolean
    if to print the lambda value
  """
  print ""
  if print_lam : print "%6s"%lam,
  print " %12s"%ene,
  if print_uncert : print " %8s"%std,
  print ""

def _parse_folder(path,res_tem,skip,maxread) :
  """ 
  Parse a number of ProtoMS result files and calculate the ensemble average of the gradient
  
  Parameters
  ----------
  path : string 
    the directory where the result files are located
  res_tem : string
    the prefix of the result files
  skip : int
    number of snapshots to skip
  maxread : int
    maximum number of snapshots to read

  Returns
  -------
  float
    the lambda value
  float
    the average gradient
  float
    the standard deviation of the average gradient
  """
   
  # List all results files and sort them
  filenames = glob.glob(path+"/%s*"%res_tem)
  if len(filenames) > 1 : filenames.sort()

  gradsum = 0.0
  gradsum2 = 0.0
  n = 0
  lam = None

  for f in filenames:
    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=f,skip=skip,readmax=maxread)
    for snapshot in results_file.snapshots:
      # what you want here is to get the current lambda and the gradient for this snapshot, so that you can do (9) 
      if hasattr(snapshot,"gradient") :
        gradient = snapshot.gradient
        lam = snapshot.lam
      else :
        lam = snapshot.lam
        lamB = snapshot.lamb
        lamF = snapshot.lamf
        dGB = snapshot.backfe
        dGF = snapshot.forwfe
        # Calculate and return the gradient
        deltalam = max(lam-lamB,lamF-lam)
        # This is needed for the end-points
        if (lam < 0.0001):
          dGB = -dGF
        if (lam > 0.9999):
          dGF = -dGB
        gradient = (dGF - dGB) / (2*deltalam)      
      if gradient != None :
        gradsum = gradsum + gradient
        gradsum2 = gradsum2 + gradient*gradient
        n = n + 1
    # Computes ensemble average and standard error
    av = gradsum / n
    std = 0.0
    if n > 2 :
      std = np.sqrt((gradsum2-av*gradsum)/(n-1)/n)
    return lam,av,std

def _calc_gradients(path,res_tem,skip,maxread,print_grad,print_uncert,print_lam) :
  """
  Calculate gradients for a number of folders
  
  Parameters
  ----------
  path : string
    the directory to find lambda folders
  res_tem : string
    the prefix of the results files
  skip : int
    number of snapshots to ignore
  maxread : int
    maximum number of snapshots to read
  print_grad : boolean
    indicate if gradients should be printed
  print_uncert : boolean
    indicate if uncertainties should be printed
  print_lam : boolean
    indicate if lambda values should be printed

  Returns
  -------
  numpy array
    all lambda values
  numpy array
    gradient for each lambda value
  numpy array
    standard deviation for each lambda value
  """

  # List all lambda folders and sort them
  paths = glob.glob(path+"/lam-*")
  paths.sort()

  # Process all lambda folders
  gradients = []
  stds = []
  lambdas = []
  for path in paths :
    (lam,grad,std) = _parse_folder(path,res_tem,skip,maxread)
    if print_grad : _print_ene(lam,grad,std,print_uncert,print_lam)
    gradients.append(grad)
    lambdas.append(lam)
    stds.append(std)
  return np.array(lambdas),np.array(gradients),np.array(stds)

def ti(path,res_tem,skip,maxread,print_grad,print_pmf,print_uncert,print_lam) :
  """
  Do thermodynamic integration
  
  Parameters
  ----------
  path : string 
    directory where all lambda folders are
  res_tem : string 
    result files
  skip :  int
    number of snapshots to skip
  maxread : int 
    maximum number of snapshots to read
  print_grad : boolean
    if the gradient should be printed to standard output
  print_pmf : boolean 
    if the PMF should be calculate and printed to standard output
  print_uncert : boolean
    if to print the uncertainty
  print_lam : boolean
    if to print the lambda value

  Returns
  -------
  numpy array
    all lambda values
  numpy array
    gradient for each lambda value
  numpy array
    standard deviation of the gradient for each lambda value
  numpy array
    pmf at each lambda value
  numpy array
    standard deviation of the pmf at each lambda value
  """

  # Calculate the gradient
  if print_grad :
    _print_head("lambda","gradient","std",print_uncert,print_lam)
  lambdas,gradients,stds = _calc_gradients(path,res_tem,skip,maxread,print_grad,print_uncert,print_lam)

  
  # Calculate and print the PMF 
  pmf = np.zeros(gradients.shape)
  pmf_std = np.zeros(gradients.shape)
  if print_pmf : _print_head("lambda","PMF","std",print_uncert,print_lam)
  pmf[0] = 0.0
  w = 0.5*(lambdas[0]+lambdas[1])
  pmf_std[0] = w**2*stds[0]**2
  if print_pmf : _print_ene(lambdas[0],pmf[0],np.sqrt(pmf_std[0]),print_uncert,print_lam)
  # Trapezium integration
  for i in range(1,len(lambdas)) :
    h = lambdas[i]-lambdas[i-1]
    pmf[i] = pmf[i-1] + h*(gradients[i]+gradients[i-1])/2.0
    if i == len(lambdas) - 1 :
      w = 1.0 - 0.5*(lambdas[i]+lambdas[i-1])
    else :
      w = 0.5*(lambdas[i+1]-lambdas[i-1])      
    pmf_std[i] = pmf_std[i-1] + w**2*stds[i]**2
    if print_pmf : _print_ene(lambdas[i],pmf[i],np.sqrt(pmf_std[i]),print_uncert,print_lam)
  
  return lambdas,gradients,stds,pmf,np.sqrt(pmf_std)
#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate free energy from thermodynamic integration")
  parser.add_argument('-d','--directory',help="the root directory that contains all the output files of the simulation. Default is cwd.",default="./")
  parser.add_argument('-r','--results',help="the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
  parser.add_argument('-pg','--print-grad',dest='printGrad',action='store_false',help="turns off printing of gradient",default=True)
  parser.add_argument('-pp','--print-pmf',dest='printPMF',action='store_false',help="turns off printing of PMF",default=True)
  parser.add_argument('-pu','--print-uncert',dest='printUncert',action='store_false',help="turns off printing of uncertainties",default=True)
  parser.add_argument('-pl','--print-lam',dest='printLam',action='store_false',help="turns off printing of lambda-values",default=True)
  parser.add_argument('-gr','--plot-grads',dest='plotGrad',action='store_false',help="turns off producing plot of gradients",default=True)
  args = parser.parse_args()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1
  # Do thermodynamic integration
  lambdas,gradients,grad_std,pmf,pmf_std = ti(args.directory,args.results,args.skip,args.max,args.printGrad,
                                               args.printPMF,args.printUncert,args.printLam)

  # Plot the gradient
  if args.plotGrad :
    pl.plot(lambdas,gradients,color='c')
    pl.xlabel('$\lambda$')
    pl.ylabel('$\Delta G / \Delta \lambda$ (kcal)')
    pl.savefig('gradients.png')
   
 
