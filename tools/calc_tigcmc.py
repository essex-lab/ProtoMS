# Author: Gregory Ross


"""
Program to do thermodynamic integration

This module defines a single public function:
tigcmc

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

from scipy.stats import spearmanr

import simulationobjects
import calc_gci
import calc_ti

logger = logging.getLogger('protoms')


def _calc_gradmatrix(path,res_tem,skip,maxread,verbose,numkind,useanalytical):
    """ 
    Parse a number of ProtoMS result files and calculate the ensemble average of the gradient for each
    lambda Adams value pair. 
  
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
    # Find lambdas and sort them
    paths = glob.glob(os.path.join(path,"lam-*"))
    paths.sort()
    gcfldr = glob.glob(paths[1]+"/b_*")
    
    # These matrices store the gradients, average number of waters and B values from each lambda-B value results file.
    # Lambdas will be along the columns, B values on the rows.
    gradmat = np.zeros((len(paths),len(gcfldr)))
    watmat = np.zeros((len(paths),len(gcfldr)))
    adamsmat = np.zeros((len(paths),len(gcfldr)))
    lambmat = np.zeros((len(paths),len(gcfldr)))
    stdmat = np.zeros((len(paths),len(gcfldr)))

    i = 0    # The row index of gradmat 
    j = 0    # The column index of gradmat
    # Process get gradients in lambda folders, as well as their corresponding gradients and water occupancies.
    for lamfldr in paths:
        N = []
        B = []
        gcmcfolders = glob.glob(lamfldr+"/b_*")
        gcmcfolders.sort(reverse=True)
        for gcfldr in gcmcfolders:
            (lambmat[i,j],gradmat[i,j],stdmat[i,j], watmat[i,j],adamsmat[i,j]) = calc_ti.parse_folder(gcfldr,res_tem,skip,maxread,numkind,useanalytical)
            j = j + 1        # Updating so the gradient at the next B value can be stored
        j = 0                # Restarts the column index
        i = i + 1            # Moves to the next row 
    # THE MATRICES SHOULD BE ORDERED PROPERLY!!!
    # May be okay, as final collapsed matrices MAY be properly ordered.
    return lambmat,gradmat,stdmat, watmat,adamsmat

def _calc_gcweights(lambmat,gradmat,stdmat,watmat,adamsmat,singlepath=False,hyd_nrg =-6.2,T=298.15):
 
    kT = T*0.0019872
    
    gci_energy_mat = np.zeros(np.shape(adamsmat))
    gci_weight_mat = np.zeros(np.shape(adamsmat))
    
    # Check to make sure all B values are the same
    
    # MOVE THESE TO THE READ DATA FUNCTION
    if adamsmat.std(axis=0).max() != 0.0:
        print "Error: B values between lambda folders don't match!"
        return
    if lambmat.std(axis=1).max() != 0.0:
        print "Error: Unequal lambda values found between replicas!"
        return
    
    # Looping over each lambda value and calculating the GCI free energy
    for i in range(0,np.shape(watmat)[0]):
        #print watmat[i,]
        
        # Smoothing data by fitting an artifial neural network for a given lambda. Note that the fitting is independent of the other lambdas.
        single_model, models = calc_gci.fit_ensemble(x=adamsmat[i,],y=watmat[i,],size=size,verbose=False,randstarts=10)
        
        # Based on the (best) model, what is the predicted number of waters at each B]
        smoothed_watnums = single_model.predicted

        # Calculating ideal gas transfer free energies and weights for each B value
        gc_bind_nrgs = []
        for num in smoothed_watnums:
            N_range = np.array((smoothed_watnums.min(),num))
            trans_nrg = calc_gci.insertion_pmf(N_range,single_model)[1]
            gc_bind_nrgs.append(trans_nrg - hyd_nrg*num)
            
        gci_energy_mat[i,] = np.array(gc_bind_nrgs)
        
        if singlepath == False:
            weights = np.exp(-gci_energy_mat[i,]/kT)
            gci_weight_mat[i,] =  weights/np.sum(weights)
        else:
            weights = (gci_energy_mat[i,] == gci_energy_mat[i,].min())*1.0   # Finds the B value with lowest binding free energy
            gci_weight_mat[i,] =  weights/np.sum(weights)                    # In case the minimum is not unique, give equal weight
        
    return gci_energy_mat, gci_weight_mat


def _collapse_matrices(gradmat,weightmat,lambmat,stdmat):
    gradients = (weightmat*gradmat).sum(axis=1)    
    std = np.sqrt((weightmat*(stdmat**2)).sum(axis=1))    # Weighted sum of variance. May or may not be the correct formular
    return lambmat[:,1],gradients,std


def tigcmc(path,res_tem,skip,maxread,verbose,numkind,useanalytical,singlpath) :
  """
  Do thermodynamic integration using data from GCMC-lambda replica exchange. Based on calc_ti.ti()
  
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
    calc_ti.print_head("lambda","gradient","std",verbose["uncert"],verbose["lambda"])

  (lambmat,gradmat,stdmat,watmat,adamsmat) = _calc_gradmatrix(path,res_tem,skip,maxread,verbose,numkind,useanalytical)
  (energy_mat,weightmat) = _calc_gcweights(lambmat,gradmat,stdmat, watmat,adamsmat,singlepath)
  lambdas,gradients,stds = _collapse_matrices(gradmat,weightmat,lambmat,stdmat)
  
  # Calculate and print the PMF 
  pmf = np.zeros(gradients.shape)
  pmf_std = np.zeros(gradients.shape)
  if verbose["pmf"] : calc_ti.print_head("lambda","PMF","std",verbose["uncert"],verbose["lambda"])
  pmf[0] = 0.0
  w = 0.5*(lambdas[0]+lambdas[1])
  pmf_std[0] = w**2*stds[0]**2
  if verbose["pmf"] : calc_ti.print_ene(lambdas[0],pmf[0],np.sqrt(pmf_std[0]),verbose["uncert"],verbose["lambda"])
  # Trapezium integration
  for i in range(1,len(lambdas)) :
    h = lambdas[i]-lambdas[i-1]
    pmf[i] = pmf[i-1] + h*(gradients[i]+gradients[i-1])/2.0
    if i == len(lambdas) - 1 :
      w = 1.0 - 0.5*(lambdas[i]+lambdas[i-1])
    else :
      w = 0.5*(lambdas[i+1]-lambdas[i-1])      
    pmf_std[i] = pmf_std[i-1] + w**2*stds[i]**2
    if verbose["pmf"] : calc_ti.print_ene(lambdas[i],pmf[i],np.sqrt(pmf_std[i]),verbose["uncert"],verbose["lambda"])

  if not verbose["pmf"] and verbose["total"] : 
    calc_ti.print_ene(lambdas[-1],pmf[-1],np.sqrt(pmf_std[-1]),verbose["uncert"],verbose["lambda"])  

  return lambdas,gradients,stds,pmf,np.sqrt(pmf_std)
#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate free energy from thermodynamic integration")
  parser.add_argument('-d','--directory',help="the root directory that contains all the output files of the simulation. Default is cwd.",default="./")
  parser.add_argument('-r','--results',help="the name of the file to analyse. Default is results. ",default="results_inst")
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
  parser.add_argument('--singlepath',action='store_true',default=False, help="whether to integrate over the lowest GCMC free energy B value when doing TI, instead of performing a weighted average")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1
  # Do thermodynamic integration
  verbose = {"total":True,"gradient":args.printGrad,"pmf":args.printPMF,"uncert":args.printUncert,"lambda":args.printLam}
  lambdas,gradients,grad_std,pmf,pmf_std= tigcmc(args.directory,args.results,args.skip,args.max,verbose,args.numerical,args.analytical,args.singlepath)

  # Do the fit
  if args.fitPMF :
    fitcoef = fit_pmf(lambdas,pmf)

  # Plot the gradient
  if args.plotGrad :
    pl.plot(lambdas,gradients,color='c')
    pl.xlabel('$\lambda$')
    pl.ylabel('$\Delta G / \Delta \lambda$ (kcal)')
    pl.savefig('gradients.png')
   
 
