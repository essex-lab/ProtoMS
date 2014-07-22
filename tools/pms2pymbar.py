# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to prepare ProtoMS output for pymbar
"""

import glob
import os

import numpy as np

import simulationobjects

def _parse_folder(path,res_tem,skip,maxread) :
  """ 
  Parse a number of ProtoMS result files and extract MBAR energies
  
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
   the current lambda value
  numpy array
    the total energy at each lambda value for each snaphot 
  """
   
  filename = os.path.join(path,res_tem)
  results_file = simulationobjects.ResultsFile()
  results_file.read(filename=filename,skip=skip,readmax=maxread) 
  nlambda = len(results_file.snapshots[0].feenergies)
  totalenergies = np.zeros([len(results_file.snapshots),nlambda])

  for i,snapshot in enumerate(results_file.snapshots) :
    for j,lam in enumerate(sorted(snapshot.feenergies.keys())) :
      totalenergies[i,j] = snapshot.feenergies[lam]
  return results_file.snapshots[0].lam,totalenergies
    
def _extract_energies(path,res_tem,skip,maxread) :
  """
  Extract total energies from a number of folders
  
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

  Returns
  -------
  list of float 
    the lambda values
  list of numpy array
    all total energies
  list of strings
    the path of the lambda folders
  """

  # List all lambda folders and sort them
  paths = glob.glob(path+"/lam-*")
  paths.sort()

  total_energies = []
  lambdas = []
  for path in paths :
    lam,ene = _parse_folder(path,res_tem,skip,maxread)
    lambdas.append(lam)
    total_energies.append(ene)
    
  return lambdas,total_energies,paths

def _write_file(filename,fromname,lam,energies) :
  """
  Write the total energies to a file that can be parsed by pymar
  
  Parameters
  ----------
  filename : string
    the name of the file
  fromname : string
    the original name of the result file
  lam : float
    the current lambda
  energies : numpy array
    energies at all lambda values and snapshots
  """
  
  with open(filename,"w") as f :
    f.write("# ProtoMS FE data extracted from %s\n"%fromname)
    f.write("#\n")
    f.write("# Energies in kcal/mol\n")
#    f.write("# Found timestep = 1.0 ps\n")
    f.write("# Found lambda = %.3f\n"%lam)
#    f.write("# Found bar_intervall = 1\n")
    f.write("# Found total of %d l-values\n"%energies.shape[1])
    f.write("# Energy values:\n")
    for i,ene in enumerate(energies) :
      f.write("%12.1f 0.0000 %s 1.000\n"%(i," ".join("%20.8E"%e for e in ene)))

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to extract ProtoMS results for pymbar")
  parser.add_argument('-d','--directory',help="the root directory that contains all the output files of the simulation. Default is cwd.",default="./")
  parser.add_argument('-r','--results',help="the name of the file to analyse. Default is results. ",default="results")
  parser.add_argument('-o','--out',help="the name of the file to write. Default is pymbar_energy. ",default="pymbar_energy")
  parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)#
  parser.add_argument('-t','--temperature',type=float,help="the simulation temperature in degrees. Default is 25.0",default=25.0)
  parser.add_argument('--run',action='store_true',help="whether to run pymbar",default=False)
  args = parser.parse_args()  

  lambdas,energies,paths = _extract_energies(args.directory,args.results,args.skip,args.max)
  for lam,ene,path in zip(lambdas,energies,paths) :
    filename = os.path.join(path,args.out)
    _write_file(filename,args.results,lam,ene)
    print "Writing pymbar energies to: %s"%filename

  if args.run :
    try :
      import pymbar
    except :
      print "Could not import pymbar! Please make sure that pymbar is in your Python path."
      quit()
    
    RT = 1.9872041*(args.temperature+273.15)/1000.00
    nstates = len(energies)
    N_k = np.zeros(nstates) 
    for i,e in enumerate(energies) : N_k[i] = e.shape[0]
    u_kln = np.zeros([nstates,nstates,N_k.max()]) 
    for i,e in enumerate(energies) :
      u_kln[i,:,:e.shape[0]] = e.T / RT
    MBAR = pymbar.MBAR(u_kln,N_k)
    (Deltaf_ij, dDeltaf_ij) = MBAR.getFreeEnergyDifferences(uncertainty_method='svd-ew')
    print "MBAR estimate: %-6.2f +- %-6.2f"%(Deltaf_ij[0,nstates-1]*RT,dDeltaf_ij[0,nstates-1]*RT)
