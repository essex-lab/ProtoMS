# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to plot the theta distribution, resulting from a JAWS
stage one simulation

This module defines the following public functions:

Can be executed from the command line as a stand-alone program
"""

import os
import logging

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pylab as plt

import simulationobjects

logger = logging.getLogger('protoms')

def solname_to_id(molecule,restart='restart') :
  """
  Search for the residue name of a solute in the
  restart file and find out the correspondent 
  solute id (or ids, if there are several copies)

  Parameters
  ----------
  molecule : string
    the residue name of the JAWS
    molecule
  restart : string, optional
    the name of the restart file
  
  Returns
  -------
  list :
    the list of solute / gcsolute ids that
    correspond to that residue name
  string :
    the molecule type. Either "solute", "gcsolute"
    or "both"
  """

  restart_obj = simulationobjects.RestartFile(filename=restart)
  sol_list = [key for key in restart_obj.solutesdic if restart_obj.solutesdic[key] in molecule.lower()]
  gcsol_list = [key for key in restart_obj.gcsolutesdic if restart_obj.gcsolutesdic[key] in molecule.lower()]

  if not sol_list :
    return gcsol_list, "gcsolute"
  elif not gcsol_list :
    return sol_list, "solute"
  else :
    return [gcsol_list,sol_list], "both"


def find_solute_theta(molecule,restart='restart',results='results') :
  """
  Given a solute / gcsolute residue name, find
  and store its corresponding values of theta through
  the simulation.

  Parameters
  ----------
  molecule : string
    the residue name of the JAWS
    molecule
  restart : string, optional
    the name of the restart file
  results : string, optional
    the name of the results file

  Returns
  -------
  dict : 
    the theta values per solute id
    in snapshot order 
  """

  id_list,mol_type = solname_to_id(molecule,restart=restart)
  results_obj = simulationobjects.ResultsFile()
  results_obj.read(filename=results)
  thetas = dict([[i,[]] for i in id_list])
  for snap in results_obj.snapshots :
    if mol_type is "solute" :
      res_thetavals = snap.thetasolvals
    elif moltype is "gcsolute" :
      res_thetavals = snap.thetavals
    for ind,res_theta in enumerate(res_thetavals) :
      thetas[ind+1].append(res_theta)
  return thetas

def plot_theta_dist(thetas) :
  """
  Plot a distribution of thetas
  """
  def set_axes() :
    plt.xlabel(r"$\theta$",fontsize=20)
    plt.ylabel("Frequency",fontsize=15)

  thetas_values = np.array(thetas.values(),float)

  n, bins, patches = plt.hist(np.reshape(thetas_values,newshape=-1),bins=100,facecolor=simulationobjects.color(0),histtype='stepfilled')
  set_axes()
  plt.savefig("theta_dist_all.png")
  plt.clf()

  n, bins, patches = plt.hist(np.transpose(thetas_values),bins=100,histtype="step")
  set_axes()
  plt.savefig("theta_dist_each.png")
  
  

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to plot the theta distribution of a given molecule, result from a JAWS simulation")
  parser.add_argument('-r','--results',help="the name of the results file. Deafult='results'",default='results')
  parser.add_argument('-s','--restart',help="the replica values to plot. Default='restart'",default='restart')
  parser.add_argument('-m','--molecule',help="the residue name of the JAWS molecule. Default='WAT'",default="WAT")
  args = parser.parse_args()

  thetas_dic = find_solute_theta(args.molecule,args.restart,args.results)

  plot_theta_dist(thetas_dic)
