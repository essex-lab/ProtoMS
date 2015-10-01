# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to plot the theta distribution, resulting from a JAWS
stage one simulation

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

from calc_clusters import cluster_coords

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
  sol_list = [key for key in restart_obj.solutesdic if restart_obj.solutesdic[key].lower() in molecule.lower()]
  gcsol_list = [key for key in restart_obj.gcsolutesdic if restart_obj.gcsolutesdic[key].lower() in molecule.lower()]

  if not sol_list :
    return gcsol_list, "gcsolute"
  elif not gcsol_list :
    return sol_list, "solute"
  else :
    return [gcsol_list,sol_list], "both"


def find_solute_theta(molecule,restart='restart',results='results',skip=0,print_average=False) :
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
  skip : int, optional
    the number of snapshots to skip

  Returns
  -------
  dict : 
    the theta values per solute id
    in snapshot order 
  """

  id_list,mol_type = solname_to_id(molecule,restart=restart)
  results_obj = simulationobjects.ResultsFile()
  results_obj.read(filename=results,skip=skip)
  thetas = dict([[i,[]] for i in id_list])
  for snap in results_obj.snapshots :
    if mol_type is "solute" :
      res_thetavals = snap.thetasolvals
    elif mol_type is "gcsolute" :
      res_thetavals = snap.thetavals
    for ind,res_theta in enumerate(res_thetavals) :
      if ind+1 in thetas.keys() : thetas[ind+1].append(res_theta)
  return thetas

def thres_and_mean(dict_th, threshold=0.95) :
  """
  Calculate the mean and the proportion of
  values above a threshold from a dictionary of values

  Parameters
  ----------
  dict_th : dict
    a dictionary with several keys
    which 'values' are lists of floats
    to calculate mean and proportion
    above threshold
  threshold: optional, float
    the threshold above which the
    proportion of values will be calculated

  Return
  ------
  float
    mean theta
  float
    fraction of thetas above threshold
    
  """

  thetas_values = np.array(dict_th.values(),float)

  mean = np.mean(thetas_values)

  above_thres = float(thetas_values[np.where(thetas_values > threshold)].size) / float(thetas_values.size)

  return  mean, above_thres


def plot_dist(dict_th,plotname="theta_dist",xlab="theta",xlimit=None,bigaxes=False) :
  """
  Plot a distribution of dictonary values
  in form of lists

  Parameters
  ----------
  dict_th : dict
    a dictionary with several keys
    which 'values' are lists of floats
    to plot in a histogram
  plotname : string, optional
    the beginning of the file name
    for the plots produced
  xlab : string, optional
    what to write in the x axes label
    of the histogram
  xlimit : list, optional
    the limits of the x-axes values.
    Set to [0.0,0.1] if 'theta' in xlab
  bigaxes : boolean, optional
    whether to set the label and axes
    of the plots to a bigger font

  Return
  ------
  Nothing
    the plots will be saved in png files
  """

  def set_axes(xname,yname,ax) :
    fontsizex = fontsizey = 15
    fontweight = None
    if bigaxes :
      fontsizex = fontsizey = 24
      xticksize = yticksize = 17
      plt.setp(ax.get_xticklabels(), fontsize=xticksize)
      plt.setp(ax.get_yticklabels(), fontsize=yticksize)
    if '$' in xname : fontsizex = fontsizex*1.5
    if '$' in yname : fontsizey = fontsizey*1.5
    plt.xlabel(xname,fontsize=fontsizex)
    plt.ylabel(yname,fontsize=fontsizey)
    xends = xlimit
    if "theta" in xname and not xlimit : xends = [0.0,1.0]
    if xends : plt.xlim(tuple(xends))

  if xlab == "theta" : xlab = r'$\theta$'
      

  thetas_values = np.array(dict_th.values(),float)

  ax = plt.figure().add_subplot(1,1,1)
  set_axes(xlab,"frecuency",ax)
  n, bins, patches = ax.hist(np.reshape(thetas_values,newshape=-1),bins=100,facecolor=simulationobjects.color(0),histtype='stepfilled')
  plt.savefig(plotname+"_all.png")
  plt.clf()

  ax = plt.figure().add_subplot(1,1,1)
  set_axes(xlab,"frecuency",ax)
  n, bins, patches = ax.hist(np.transpose(thetas_values),bins=100,histtype="step")
  plt.savefig(plotname+"_each.png")
  plt.clf()

  ax = plt.figure().add_subplot(1,1,1)
  set_axes("snapshot",xlab,ax)
  for i in dict_th :
    ax.plot(dict_th[i])
  plt.savefig(plotname+"_snapshot.png")
  
    

def extract_theta_pdb (thetas_dic, theta_range, pdbfile, outpdb, residue) :
  """
  Extract the molecules which theta is
  within a given range
  
  Parameters
  ----------
  thetas_dic : dict
    the theta values per solute id
    in snapshot order 
  theta_range : list
    the lower and upper limit of the
    theta range
  pdbfile : string
    the name of the pdb file product
    of the ProtoMS simulation
  outpdb : string
    the name of the pdbfile where
    the extracted molecules will be
    printed

  Return
  ------
  PDBFile :
    the PDBFile object with the
    molecules within the given theta range
  """

  outpdb = "%s_%s-%s.pdb"%(outpdb,theta_range[0],theta_range[1])

  try :
    theta_range = [float(limit) for limit in theta_range]
  except :
    raise simulationobjects.SetupError("The limits of the theta range could not be undestood")

  pdbin_obj = simulationobjects.PDBSet()
  pdbin_obj.read(filename=pdbfile)


  pdbout_obj = simulationobjects.PDBFile()
  pdbout_obj.header = "HEADER   %s molecules found with theta between %3.3f and %3.3f\n"%(residue, theta_range[0], theta_range[1])
  
  for snapshot, pdbin_file in enumerate(pdbin_obj.pdbs) :
    dic_inresidues = pdbin_file.residues
    if simulationobjects.is_solvent(residue) : dic_inresidues = pdbin_file.solvents 
    my_inresidues = [dic_inresidues[reskey] for reskey in dic_inresidues if dic_inresidues[reskey].name in residue]
    for ind, myres in enumerate(my_inresidues) :
      if theta_range[0] <= float(thetas_dic[thetas_dic.keys()[ind]][snapshot]) <= theta_range[1] :
        pdbout_obj.residues[(snapshot+1)*(ind+1)] = myres

  pdbout_obj.write(filename=outpdb)
          

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
  parser.add_argument('-p','--plotname',help="the start of the filename for the plots generated. Default='theta_dist'",default="theta_dist")
  parser.add_argument('-a','--average',action="store_true",help="whether to print an average of all theta values. Default = False",default=False)
  parser.add_argument('-tr','--thetarange',nargs=2,help="the range of thetas which correponding molecules will be extracted, specified by a lower and an upper limit. By default, no molecules will be extracted.")
  parser.add_argument('-pdb','--pdbfile',help="the pdb file result of protoms. Only required if when molecules are extracted for a range of thetas. Default = 'all.pdb'",default='all.pdb')
  parser.add_argument('-op','--outpdb',help="the name of the file where the molecules for a given range of thetas are to be extracted. Default = 'out_thetas'",default='out_thetas')
  parser.add_argument('-th','--threshold',help="the proportion of thetas above this threshold will be calculated, as well as the total mean theta. By default this function is not active",default=None)
  parser.add_argument('--skip',help="the number of results snapshots to skip, Default = 0",default="0")
  parser.add_argument('--bigaxes',action="store_true",help="whether to set the axes in the plots to big fontzise. Default=False",default=False)
  args = parser.parse_args()

  try :
    skip_steps = int (args.skip)
  except :
    simulationobjects.SetupError("Snapshots to skip needs to be an integer. The argument %s could not be interpreted."%args.skip)

  logger = simulationobjects.setup_logger("plot_theta.log")

  thetas_dic = find_solute_theta(args.molecule,args.restart,args.results,skip_steps,args.average)

  plot_dist(thetas_dic,plotname=args.plotname,bigaxes=args.bigaxes)

  if args.threshold :
    try :
      threshold = float(args.threshold)
    except :
      simulationobjects.SetupError("Threshold value %s could not be interpreted as float."%args.threshold)

    mean_theta, above_threshold = thres_and_mean(thetas_dic,threshold)
    print "Mean theta: %.3f"%mean_theta
    print "Proportion of values above threhold %.3f: %.2f%%"%(threshold,above_threshold*100)

  if args.thetarange :
    extract_theta_pdb(thetas_dic,args.thetarange,args.pdbfile,args.outpdb,args.molecule)









