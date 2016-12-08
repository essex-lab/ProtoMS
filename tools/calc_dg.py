# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to calculate free energies using TI, BAR and MBAR
"""

import os
import logging

import numpy as np

import simulationobjects
import calc_ti
import calc_bar
import pms2pymbar

import matplotlib.pyplot as plt

logger = logging.getLogger('protoms')

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate free energy from TI/BAR/MBAR")
  parser.add_argument('-d','--directories',nargs="+",help="the root directory that contains all the output files of the simulation. Default is cwd.",default=["./"])
  parser.add_argument('-r','--results',help="the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
  parser.add_argument('-t','--temperature',type=float,help="the simulation temperature in degrees. Default is 25.0",default=25.0)
  parser.add_argument('-b','--nboots',type=int,help="the number of bootstrap samples",default=100)
  parser.add_argument('-gr','--plot-grads',dest='plotGrad',action='store_true',help="turns on producing plot of gradients",default=False)
  parser.add_argument('-pg','--print-grad',dest='printGrad',action='store_false',help="turns off printing of gradient",default=True)
  parser.add_argument('-pe','--print-each',dest='printEach',action='store_false',help="turns off printing of each free energy",default=True)
  parser.add_argument('--analytical',action='store_true',help="turns on use of analytical gradients",default=False)
  parser.add_argument('--numerical',choices=["both","back","forw"],default="both",help="the kind of numerical gradient estimator")
  parser.add_argument('-e','--estimator',nargs="+",choices=["ti","bar","mbar"],default=["ti","bar","mbar"],help="the type of estimator to use")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1

  RT = 1.9872041*(args.temperature+273.15)/1000

  # Estimate TI for each of the output directories
  if "ti" in args.estimator :
    print ""
    verbose = {"total":False,"gradient":False,"pmf":False,"uncert":False,"lambda":False}  
    lambdas = None
    all_gradients = []
    all_grad_std = []  
    pmf = pmf_std = None
    for directory in args.directories :
      lambdas,gradients,grad_std,pmf,pmf_std = calc_ti.ti(directory,args.results,args.skip,args.max,verbose,args.numerical,args.analytical)
      all_gradients.append(gradients)
      all_grad_std.append(grad_std)
      if args.printGrad and len(args.directories) == 1 :
        print "\n%6s %12s %8s"%("lambda","gradient","std")
        for i in range(len(lambdas)) :
          print "%-6.3f %12.4f %8.4f"%(lambdas[i],gradients[i],grad_std[i])
        print ""
      if args.printEach or len(args.directories) == 1:        
        print "TI (%s) : %.3f +- %.3f"%(directory,pmf[-1],pmf_std[-1])
    if len(all_gradients) == 1 :    
      if args.plotGrad :
        plt.plot(lambdas,all_gradients[0],color='c')
    else :
      # Calculates the average gradient and do the integration over that
      gradients = np.array(all_gradients).mean(axis=0)
      grad_std = np.array(all_gradients).std(axis=0)/np.sqrt(len(all_gradients))
      pmf = np.zeros(gradients.shape)
      pmf_std = np.zeros(gradients.shape)
      w = 0.5*(lambdas[0]+lambdas[1])
      pmf_std[0] = w**2*grad_std[0]**2
      if args.printGrad : 
        print "\n%6s %12s %8s"%("lambda","gradient","std")
        print "%-6.3f %12.4f %8.4f"%(lambdas[0],gradients[0],grad_std[0])

      for i in range(1,len(lambdas)) :
        if args.printGrad : print "%-6.3f %12.4f %8.4f"%(lambdas[i],gradients[i],grad_std[i])
        h = lambdas[i]-lambdas[i-1]
        pmf[i] = pmf[i-1] + h*(gradients[i]+gradients[i-1])/2.0
        if i == len(lambdas) - 1 :
          w = 1.0 - 0.5*(lambdas[i]+lambdas[i-1])
        else :
          w = 0.5*(lambdas[i+1]-lambdas[i-1])      
        pmf_std[i] = pmf_std[i-1] + w**2*grad_std[i]**2
      print "TI : %.3f +- %.3f"%(pmf[-1],np.sqrt(pmf_std[-1]))
      if args.plotGrad :
        plt.errorbar(lambdas,gradients,yerr=grad_std,color='c')
    if args.plotGrad :
      plt.xlabel('$\lambda$')
      plt.ylabel('$\Delta G / \Delta \lambda$ (kcal)')
      plt.savefig('gradients.png')  

  # Estimate BAR for each of the output directories
  if "bar" in args.estimator :
    print ""
    verbose = {"total":False,"window":False,"uncert":False,"lambda":False}
    if len(args.directories) :
      args.nboots = 2
    all_dgs = []
    for directory in args.directories :
      lambdas,windg,std = calc_bar.bar(directory,args.results,args.skip,args.max,RT,verbose,args.nboots)
      all_dgs.append(windg.sum())
      if args.printEach or len(args.directories) == 1:
        print "BAR (%s): %.3f +- %.3f"%(directory,windg.sum(),np.sqrt(np.sum(std**2)))
    if len(all_dgs) > 1 :
      dg = np.array(all_dgs).mean()
      std = np.array(all_dgs).std() / np.sqrt(len(all_dgs))
      print "BAR: %.3f +- %.3f"%(dg,std) 

  # Estimate MBAR for each of the output directories
  if "mbar" in args.estimator :
    print ""
    try :
      import pymbar
    except :
      print "Could not import pymbar! Please make sure that pymbar is in your Python path."
      quit()

    all_dgs = []
    for directory in args.directories : 
      lambdas,energies,paths = pms2pymbar.extract_energies(directory,args.results,args.skip,args.max)
      try:
        resp = pms2pymbar.mbar(lambdas,energies,RT)
      except np.linalg.LinAlgError:
        print "MBAR (%s): nan +- nan" % directory
        continue
      (Deltaf_ij, dDeltaf_ij) = (resp[0],resp[1])
      if args.printEach or len(args.directories) == 1:
        print "MBAR (%s): %.3f +- %.3f"%(directory,Deltaf_ij[0,-1]*RT,dDeltaf_ij[0,-1]*RT)
      all_dgs.append(Deltaf_ij[0,-1]*RT)
    if len(all_dgs) > 1 :
      dg = np.array(all_dgs).mean()
      std = np.array(all_dgs).std() / np.sqrt(len(all_dgs))
      print "MBAR: %.3f +- %.3f"%(dg,std) 
