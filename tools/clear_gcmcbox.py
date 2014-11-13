# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routine to remove solvent molecules from a GCMC/JAWS-1 simulation box

This module defines a single public function:
clear_gcmcbox

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np

import simulationobjects

logger = logging.getLogger('protoms')

def clear_gcmcbox(gcmcbox,waters) :
  """
  Removes solvent molecules from the GCMC/JAWS-1 box as they will not be 
  able to move during the simulation
  
  Parameters
  ----------
  gcmcbox : string or PDBFile object
    the gcmcbox
  waters : string or PDBFile object
    the water molecule to check

  Returns
  -------
  int 
    the number of removed water molecules
  PDBFile
    the cleared solvation box
  """
  
  logger.debug("Running clear_gcmcbox with arguments: ")
  logger.debug("\tgcmcbox = %s"%gcmcbox) 
  logger.debug("\twaters = %s"%waters) 
  logger.debug("This will remove solvent molecules within the GCMC/JAWS box")

  extend = 1.0					# The amount to clear in excess of the box limits, as only the oxygen atoms of the water molecules are used to decide whether the water molecule is in the box or not.  
  if isinstance(gcmcbox,basestring) :
    gcmcbox = simulationobjects.PDBFile(filename=gcmcbox)
  if isinstance(waters,basestring) :
    waters = simulationobjects.PDBFile(filename=waters)
    
  # Try to find box information in the header  
  box = simulationobjects.find_box(gcmcbox)
  if box is None :
    # if that fails, take the extent of the PDB structure
    box = gcmcbox.getBox()  
  if "origin" not in box :
    box["origin"] = box["center"] - box["len"]/2.0
  box_max = box["origin"] + box["len"]
  box_min = box["origin"]

  # Remove waters that are inside the GCMC/JAWS-1 box
  nrem = 0
  removethese = []
  for soli in waters.solvents :
    xyz = waters.solvents[soli].atoms[0].coords
    if np.all(xyz < (box_max + extend)) and np.all(xyz > (box_min - extend)) :
      logger.debug("Removing water %d from %s"%(soli,waters))
      nrem = nrem + 1
      removethese.append(soli)
  for soli in removethese : del waters.solvents[soli]
  logger.info("Removed %d water molecules from %s that were inside the GCMC/JAWS box %s"%(nrem,waters,gcmcbox))
  
  return nrem,waters 

if __name__ == "__main__":

  import argparse
  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to remove water molecules from a GCMC/JAWS-1 box")
  parser.add_argument('-b','--box',help="the name of the PDB-file containing the box.")
  parser.add_argument('-s','--solvation',help="the name of the PDB-file containing the solvation waters")
  parser.add_argument('-o','--out',help="the name of the output PDB-file",default="cleared_box.pdb")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("clear_gcmcbox_py.log")

  if args.solvation is None : 
    print "No pdb with solvent provided. Nothing to do! Exiting."
  
  nrem,cleared_box = clear_gcmcbox(args.box,args.solvation)
  cleared_box.write(args.out)
