# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to do BAR

This module defines a single public function:
bar

Can be executed from the command line as a stand-alone program
"""

import glob
import os
import logging

import numpy as np

import simulationobjects

logger = logging.getLogger('protoms')

def _print_ene(lam0,lam1,ene,std,print_uncert,print_lam) :
  """ 
  Print energy row to standard output
  
  Parameters
  ----------
  lam0 : float
    the first lambda value
  lam1 : float
    the second lambda value
  ene : float
    the energy value
  std : float
    the uncertainty
  print_uncert : boolean
    if to print the uncertainty
  print_lam : boolean
    if to print the lambda value
  """
  if print_lam and lam0 is not None : 
    print "%-6.2f%-6.2f"%(lam0,lam1),
  elif lam0 is None :
    print "%12s"%"Total:",
  print " %12.4f"%ene,
  if print_uncert : print " %8.4f"%std,
  print ""

def _print_head(lam0,lam1,ene,std,print_uncert,print_lam) :
  """ 
  Print energy header to standard output

  Parameters
  ----------
  lam0 : string
    the header of the first lambda value
  lam1 : string
    the header of the second lambda value
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
  if print_lam : print "%-6s%-6s"%(lam0,lam1),
  print " %12s"%ene,
  if print_uncert : print " %8s"%std,
  print ""

def _parse_folder(path,res_tem,skip,maxread) :   
  """ 
  Parse a number of ProtoMS result files and extract work values

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
  numpy array
    the backward work values
  numpy array
    the forward work values
  """

  # List all results files and sort them
  filenames = glob.glob(path+"/%s*"%res_tem)
  if len(filenames) == 1 : # This indicates that we deal with ProtoMS3.0 results files

    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=filenames[0],skip=skip,readmax=maxread)    
    workF = []
    workB = []
    for snapshot in results_file.snapshots :
      lam = snapshot.lam
      lamvals = sorted(snapshot.feenergies.keys())
      ene_thislam = snapshot.feenergies[lam]
      lamidx = lamvals.index(lam)
      if lamidx == 0 :
        dB = 0.0
      else :
        dB = ene_thislam - snapshot.feenergies[lamvals[lamidx-1]]
      if lamidx == len(lamvals)-1 :
        dF = 0.0
      else :
        dF = ene_thislam - snapshot.feenergies[lamvals[lamidx+1]]
      workF.append(dF)
      workB.append(dB)

  else :
    filenames.sort()
    # Select some of the snapshots
    filenames = filenames[skip+1:min(maxread+skip+1,len(filenames))]

    # Process each result files
    workF = []
    workB = []
    for filename in filenames :
      for line in open(filename,'r').readlines() :
        if line.find("dG") == 0 :
          cols = line.strip().split()
          lam = float(cols[1])
          dF = float(cols[7])-float(cols[11])
          dB = float(cols[7])-float(cols[9])
          workF.append(dF)
          workB.append(dB)
          break
  
  return lam,np.array(workB),np.array(workF)

def _do_bar(work_back,work_forw,RT,nboots) :

  """
  Calculate BAR for an individual window

  Parameters
  ----------
  work_back : numpy array 
    the backwards/reverse work values
  work_forw : numpy array
    the forward work values
  RT : float 
    the thermal temperature
  nboots : int
    the number of bootstrap samples

  Returns
  -------
  float
    the free energy difference
  float
    the uncertainty of the free energy difference
  """

  # Fermi function
  def fermi(x) : 
    return 1/(1+np.exp(x))

  # Iterative BAR estimate
  def iter_bar(workB,workF,dG,RT) :

    summa = 0.0
    beta=1/RT
    NF = workF.shape[0]
    NR = workB.shape[0]
    NTot = NF + NR
    M = np.log(float(NF)/float(NR))
    C = dG + M*RT
    for i in range(50) :
      summa = 0
      for i in range(NF) :
        summa = summa + fermi(beta*(workF[i]+C))
      lFF = np.log(summa)
      summa = 0
      for i in range(NR) :
        summa = summa + fermi(beta*(workB[i]-C))
      lFR = np.log(summa)
      dG = (lFF-lFR-M+beta*C) * RT
      C = dG + M*RT
  
    return dG

  # Calculate the free energy iteratively
  dg = iter_bar(work_back,work_forw,0.0,RT)

  # Bootstrap standard deviation
  dgboots = np.zeros(nboots) 
  for i in range(nboots) :
    idxb = np.random.randint(work_back.shape[0],size=work_back.shape[0])
    idxf = np.random.randint(work_forw.shape[0],size=work_forw.shape[0])
    dgboots[i] = iter_bar(work_back[idxb],work_forw[idxf],0.0,RT)
  std = dgboots.std()

  return dg,std

def _calc_windows(path,res_tem,skip,maxread,RT,verbose,nboots) :
  """
  Calculate free energies for individual windows

  Parameters
  ----------
  path : string 
    the directory to look for lambda folders
  res_tem : string 
    the prefix of results files
  skip : int
    number of snapshots to ignore
  maxread : int 
    maximum number of snapshots to process
  RT : float 
    the thermal temperature
  verbose : dictionary of booleans
    tune the printing to standard output
  nboots : int
    the number of bootstrap samples

  Returns
  -------
  numpy array
    the lambda values
  numpy array
    the free energy for each lambda
  numpy array
    the standard error of the free energy for each lambda
  """

  # List all lambda folders and sort them
  paths = glob.glob(path+"/lam-*")
  paths.sort()

  # Process all lambda folders  
  work_back = []
  work_forw = []
  lambdas = []
  for path in paths :
      lam,wb,wf = _parse_folder(path,res_tem,skip,maxread)
      work_back.append(wb)
      work_forw.append(wf)
      lambdas.append(lam)

  dgs = []
  stds = []
  for i in range(len(paths)-1) :
    dg,std = _do_bar(work_back[i+1],work_forw[i],RT,nboots)
    if verbose["window"] : _print_ene(lambdas[i],lambdas[i+1],dg,std,verbose["uncert"],verbose["lambda"])
    dgs.append(dg)
    stds.append(std)

  return np.array(lambdas),np.array(dgs),np.array(stds)

def bar(path,res_tem,skip,maxread,RT,verbose,nboots=10) :
  """
  Do Bennet Acceptance Ratio for a set of ProtoMS output files

  Parameters
  ----------
  path : string 
    directory where all lambda folders are
  res_tem : string 
    prefix of result files
  skip : int 
    number of snapshots to skip
  maxread : int 
    maximum number of snapshots to read
  RT : float 
    the thermal temperature
  verbose : dictionary of booleans
    tune the output to standard output
  nboots : int, optional
    the number of bootstrap samples

  Returns
  -------
  numpy array
    the lambda values
  numpy array
    the free energy for each lambda
  numpy array
    the standard error of the free energy for each lambda
  """

  # Calculate the gradient
  if verbose["window"] :
    _print_head("lamB","lamF","dG","std",verbose["uncert"],verbose["lambda"])
  lambdas,windg,stds = _calc_windows(path,res_tem,skip,maxread,RT,verbose,nboots)
  
  # Calculate and print the total free energy 
  if verbose["total"] :
    sum_dg = windg.sum()
    sum_std = np.sqrt(np.sum(stds**2))
    _print_ene(None,None,sum_dg,sum_std,verbose["uncert"],verbose["lambda"])

  return lambdas,windg,stds

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse
 
  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to calculate free energy from Bennett Acceptance Ratio")
  parser.add_argument('-d','--directory',help="the root directory that contains all the output files of the simulation. Default is cwd.",default="./")
  parser.add_argument('-r','--results',help="the beginning of the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
  parser.add_argument('-t','--temperature',type=float,help="the simulation temperature in degrees. Default is 25.0",default=25.0)
  parser.add_argument('-pw','--print-win',dest='printWin',action='store_false',help="turns off printing of individual windows",default=True)
  parser.add_argument('-pu','--print-uncert',dest='printUncert',action='store_false',help="turns off printing of uncertainties",default=True)
  parser.add_argument('-pl','--print-lam',dest='printLam',action='store_false',help="turns off printing of lambda-values",default=True)
  parser.add_argument('-b','--nboots',type=int,help="the number of bootstrap samples",default=100)
  args = parser.parse_args()
  
  # Setup the logger
  logger = simulationobjects.setup_logger()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1
  RT = 1.9872041*(args.temperature+273.15)/1000
  # Do BAR
  verbose = {"total":True,"window":args.printWin,"uncert":args.printUncert,"lambda":args.printLam}
  bar(args.directory,args.results,args.skip,args.max,RT,verbose,args.nboots)
   
 
