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
import logging

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pyplot as pl
import numpy as np
import calc_series

from scipy.stats import spearmanr

import simulationobjects

logger = logging.getLogger('protoms')

def print_ene(lam,ene,std,print_uncert,print_lam) :
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

def print_head(lam,ene,std,print_uncert,print_lam) :
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

def _parse_folder(path,res_tem,skip,maxread,numkind,useanalytical,auto=False) :
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
  numkind : string
    the kind of numerical gradient, should be either both, forw, or back
  useanalytical : boolean
    if to use analytical gradients
    
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
  filenames = glob.glob(os.path.join(path,"%s*"%res_tem))
  if len(filenames) > 1 : filenames.sort()
  
  gradients = []
  lam = None

  for f in filenames:
    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=f,skip=skip,readmax=maxread)
    for snapshot in results_file.snapshots:
      # what you want here is to get the current lambda and the gradient for this snapshot, so that you can do (9) 
      if useanalytical and hasattr(snapshot,"agradient") :
        gradient = snapshot.agradient
        lam = snapshot.lam
      elif numkind == "both" and hasattr(snapshot,"gradient") :
        gradient = snapshot.gradient
        lam = snapshot.lam
      else :
        if not (hasattr(snapshot,"backfe") and hasattr(snapshot,"forwfe")) :
          raise simulationobjects.SetupError("This results file does not contain any data on which one can estimate gradients.")
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
        if numkind == "back" :
          gradient = -dGB / deltalam
        elif numkind == "forw" :
          gradient = dGF / deltalam
        else :
          gradient = (dGF - dGB) / (2*deltalam)      
      if gradient != None :
        gradients.append ( gradient )

    if auto:
      gradients = gradients[calc_series.find_equilibration(gradients,np.arange(1,len(gradients)+1)):]

    # Computes ensemble average and standard error
    return lam,np.mean(gradients),np.std(gradients,ddof=1)/len(gradients)**0.5

def _calc_gradients(path,res_tem,skip,maxread,verbose,numkind,useanalytical,subd,auto=False) :
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
  verbose : dictionary of booleans
    tune the printing to standard output
  numkind : string
    the kind of numerical gradient, should be either both, forw, or back
  useanalytical : boolean
    if to use analytical gradients
  subd : string
    optional subdirectory to check for each lambda
    
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
  paths = glob.glob(os.path.join(path,"lam-*",subd))
  paths.sort()

  # Process all lambda folders
  gradients = []
  stds = []
  lambdas = []
  for path in paths :
    (lam,grad,std) = _parse_folder(path,res_tem,skip,maxread,numkind,useanalytical,auto)
    if verbose["gradient"] : print_ene(lam,grad,std,verbose["uncert"],verbose["lambda"])
    gradients.append(grad)
    lambdas.append(lam)
    stds.append(std)
  return np.array(lambdas),np.array(gradients),np.array(stds)

def fit_pmf(lambdas,pmf,orderfit=4,upperfit=5,plotfile="fit.png"):
  """
  Fit the PMF to a polynomial function
  Parameters
  ----------
  lambdas : list
    list of lambda values
  pmf : list
    list of pmf values for each lambda
  orderfit : int, optional
    the order of the polynomial to fit your pmf to
  upperfit : int, optional
    the maximum order of the polynomial to fit to
  plotfile : string, optional
    the name of the file where to plot of the fit

  Returns
  -------
  numpy array
    the coefictients of the fit
  """

  pl.plot(lambdas,pmf,label="PMF",linewidth=2)
  pl.savefig("initial.png")

  for order in range(orderfit,upperfit) :
    fit = np.poly1d(np.polyfit(lambdas,pmf,order))
    pmfpred = np.array([fit(xval) for xval in lambdas])
    pl.plot(lambdas,pmfpred)
    ydiff = np.absolute(pmfpred - pmf)
    r = np.corrcoef(pmf,pmfpred)
    r2 = np.square(r)
    rho,p = spearmanr(pmf,pmfpred)
    if p < 0.05 :
      logger.info("\nPMF fitted to a %d order polynomial:"%order)
      logger.info(np.poly1d(fit))
      logger.info("With an R^2 of %.3f and Spearman's rho of %.3f"%(r2.item(1,0),rho))
      pl.savefig(plotfile)
      logger.info("Open %s to visualize the fit.\n"%plotfile)
      return fit.c
  logger.info("No polynomial of orders between %d and %d could be successfully fitted to your PMF" %(orderfit,upperfit))
  return None
    

def ti(path,res_tem,skip,maxread,verbose,numkind,useanalytical,subd='',auto=False) :
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
  verbose : dictionary of booleans
    tune the printing to standard output
  numkind : string
    the kind of numerical gradient, should be either both, forw, or back
  useanalytical : boolean
    if to use analytical gradients
  subd : string
    optional subdirectory to check for each lambda value

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
  if verbose["gradient"] :
    print_head("lambda","gradient","std",verbose["uncert"],verbose["lambda"])
  lambdas,gradients,stds = _calc_gradients(path,res_tem,skip,maxread,verbose,numkind,useanalytical,subd,auto)

  
  # Calculate and print the PMF 
  pmf = np.zeros(gradients.shape)
  pmf_std = np.zeros(gradients.shape)
  if verbose["pmf"] : print_head("lambda","PMF","std",verbose["uncert"],verbose["lambda"])
  pmf[0] = 0.0
  w = 0.5*(lambdas[0]+lambdas[1])
  pmf_std[0] = w**2*stds[0]**2
  if verbose["pmf"] : print_ene(lambdas[0],pmf[0],np.sqrt(pmf_std[0]),verbose["uncert"],verbose["lambda"])
  # Trapezium integration
  for i in range(1,len(lambdas)) :
    h = lambdas[i]-lambdas[i-1]
    pmf[i] = pmf[i-1] + h*(gradients[i]+gradients[i-1])/2.0
    if i == len(lambdas) - 1 :
      w = 1.0 - 0.5*(lambdas[i]+lambdas[i-1])
    else :
      w = 0.5*(lambdas[i+1]-lambdas[i-1])      
    pmf_std[i] = pmf_std[i-1] + w**2*stds[i]**2
    if verbose["pmf"] : print_ene(lambdas[i],pmf[i],np.sqrt(pmf_std[i]),verbose["uncert"],verbose["lambda"])

  if not verbose["pmf"] and verbose["total"] : 
    print_ene(lambdas[-1],pmf[-1],np.sqrt(pmf_std[-1]),verbose["uncert"],verbose["lambda"])  

  return lambdas,gradients,stds,pmf,np.sqrt(pmf_std)

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate free energy from thermodynamic integration")
  parser.add_argument('-d','--directory',help="the root directory that contains all the output files of the simulation. Default is cwd.",default="./")
  parser.add_argument('--subdir',help='optional subdirectory to check for each lamda value',default='')
  parser.add_argument('-r','--results',help="the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
  parser.add_argument('-pg','--print-grad',dest='printGrad',action='store_false',help="turns off printing of gradient",default=True)
  parser.add_argument('-pp','--print-pmf',dest='printPMF',action='store_false',help="turns off printing of PMF",default=True)
  parser.add_argument('-pu','--print-uncert',dest='printUncert',action='store_false',help="turns off printing of uncertainties",default=True)
  parser.add_argument('-pl','--print-lam',dest='printLam',action='store_false',help="turns off printing of lambda-values",default=True)
  parser.add_argument('-gr','--plot-grads',dest='plotGrad',action='store_true',help="turns on producing plot of gradients",default=False)
  parser.add_argument('-pf','--print-fit',dest='fitPMF',action='store_true',help="turns on fitting the pmf to a polynomial",default=False)
  parser.add_argument('--analytical',action='store_true',help="turns on use of analytical gradients",default=False)
  parser.add_argument('--numerical',choices=["both","back","forw"],default="both",help="the kind of numerical gradient estimator")
  parser.add_argument('--autoeqb',dest='autoeqb',action='store_true',help="use automatic equilibration detection to determine how much data is included in free energy difference")
  return parser
#
# If this is run from the command-line
#
if __name__ == '__main__' :

  args = get_arg_parser().parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1
  # Do thermodynamic integration
  verbose = {"total":True,"gradient":args.printGrad,"pmf":args.printPMF,"uncert":args.printUncert,"lambda":args.printLam}
  lambdas,gradients,grad_std,pmf,pmf_std= ti(args.directory,args.results,args.skip,args.max,verbose,args.numerical,args.analytical,args.subdir,args.autoeqb)

  # Do the fit
  if args.fitPMF :
    fitcoef = fit_pmf(lambdas,pmf)

  # Plot the gradient
  if args.plotGrad :
    pl.plot(lambdas,gradients,color='c')
    pl.xlabel('$\lambda$')
    pl.ylabel('$\Delta G / \Delta \lambda$ (kcal)')
    pl.savefig('gradients.png')
   
 
