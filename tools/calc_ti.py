# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

import glob
import os

import matplotlib.pyplot as pl
import numpy as np

import simulationobjects

#
# Calculate gradients for a number of folders
# Input:
#   path - the directory to look for lambda folders
#   res_tem - the prefix of results files
#   skip - number of snapshots to ignore
#   maxread - maximum number of snapshots to process
#   printGrad - if to print gradients to standard output
#   printUncert - if to print uncertainties to standard output
#   printLam - if to to print the lambda values to standard output
# Output:
#   an array of lambda values, an array of gradients, an array of standard errors
#
def calc_gradients(path,res_tem,skip,maxread,printGrad,printUncert,printLam) :

  # List all lambda folders and sort them
  paths = glob.glob(path+"/lam-*")
  paths.sort()

  # Process all lambda folders
  gradients = []
  stds = []
  lambdas = []
  for path in paths :
    (lam,grad,std) = parse_folder(path,res_tem,skip,maxread)
    if printGrad : printEne(lam,grad,std,printUncert,printLam)
    gradients.append(grad)
    lambdas.append(lam)
    stds.append(std)
  return (np.array(lambdas),np.array(gradients),np.array(stds))

# 
# Parse a number of ProtoMS result files and calculate the ensemble average
# of the gradient
# Input:
#   path - the directory where the result files are located
#   res_tem - the prefix of the result files
#   skip - number of snapshots to skip
#   maxread - maximum number of snapshots to read
# Output:
#   the lambda value, the gradient average and the standard error of the gradient
#
def parse_folder(path,res_tem,skip,maxread) :
   
  # List all results files and sort them
  filenames = glob.glob(path+"/%s*"%res_tem)
  filenames.sort()

  gradsum = 0.0
  gradsum2 = 0.0
  n = 0
  lam = None

  for f in filenames:
    results_file = simulationobjects.ResultsFile()
    # all_results is a list of snapshots-results
    results_file.read(filename=f,skip=skip,readmax=maxread)
    for snapshot in results_file.snapshots:
      # what you want here is to get the current lambda and the gradient for this snapshot, so that you can do (9) 
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
    return (lam,av,std)
      

#
# Print energy row to standard output
# Input:
#   lam - the lambda value
#   ene - the energy value
#   std - the uncertainty
#   printUncert - if to print the uncertainty
#   printLam - if to print the lambda value
# Output:
#   none
#
def printEne(lam,ene,std,printUncert,printLam) :

  if printLam : print "%-6.2f"%lam,
  print " %12.4f"%ene,
  if printUncert : print " %8.4f"%std,
  print ""

#
# Print energy header to standard output
# Input:
#   lamH - the header of the lambda value
#   eneH - the header of the energy value
#   stdH - the header of the uncertainty
#   printUncert - if to print the uncertainty
#   printLam - if to print the lambda value
# Output:
#   none
#
def printHead(lamH,eneH,stdH,printUncert,printLam) :

  print ""
  if printLam : print "%6s"%lamH,
  print " %12s"%eneH,
  if printUncert : print " %8s"%stdH,
  print ""

#
# Do thermodynamic integration
# Input:
#   path - directory where all lambda folders are
#   res_tem - result files
#   skip - number of snapshots to skip
#   maxread - maximum number of snapshots to read
#   printGrad - if the gradient should be printed to standard output
#   printPMF - if the PMF should be calculate and printed to standard output
#   printUncert - if to print the uncertainty
#   printLam - if to print the lambda value
#   plotGrad - if to plot the gradients
#
def ti(path,res_tem,skip,maxread,printGrad,printPMF,printUncert,printLam,plotGrad) :

  # Calculate the gradient
  if printGrad :
    printHead("lambda","gradient","std",printUncert,printLam)
  (lambdas,gradients,stds) = calc_gradients(path,res_tem,skip,maxread,printGrad,printUncert,printLam)

  # Plot the gradient
  if plotGrad :
    pl.plot(lambdas,gradients,color='c')
    pl.xlabel('$\lambda$')
    pl.ylabel('$\Delta G / \Delta \lambda$ (kcal)')
    pl.savefig('gradients.png')
  
  # Calculate and print the PMF 
  if printPMF :
    printHead("lambda","PMF","std",printUncert,printLam)
    pmf = 0.0
    w = 0.5*(lambdas[0]+lambdas[1])
    std = w**2*stds[0]**2
    printEne(lambdas[0],0.0,np.sqrt(std),printUncert,printLam)
    # Trapezium integration
    for i in range(1,len(lambdas)) :
      h = lambdas[i]-lambdas[i-1]
      pmf = pmf + h*(gradients[i]+gradients[i-1])/2.0
      if i == len(lambdas) - 1 :
        w = 1.0 - 0.5*(lambdas[i]+lambdas[i-1])
      else :
        w = 0.5*(lambdas[i+1]-lambdas[i-1])      
      std = std + w**2*stds[i]**2
      printEne(lambdas[i],pmf,np.sqrt(std),printUncert,printLam)
#
# If this is run from the command-line
#
if __name__ == '__main__' :

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
  ti(args.directory,args.results,args.skip,args.max,args.printGrad,args.printPMF,args.printUncert,args.printLam,args.plotGrad)
   
 
