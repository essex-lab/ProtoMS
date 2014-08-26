# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routine to make a GCMC/JAWS-1 simulation box from a ligand

This module defines a single public function:
make_gcmcbox

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np

import simulationobjects

logger = logging.getLogger('protoms')

def make_gcmcbox(pdb,filename,padding=2.0) :
  """
  Make a GCMC/JAWS-1 simulation box around a PDB-structure
  
  Parameters
  ----------
  pdb : PDBFile object
    the PDB structure
  filename : string
    the name of the output file
  padding : float, optional
    the amount of extra space around the ligand to add
  """
  
  logger.debug("Running make_gcmcbox with arguments: ")
  logger.debug("\tpdb      = %s"%pdb) 
  logger.debug("\tfilename = %s"%filename) 
  logger.debug("\tpadding  = %d"%padding) 
  logger.debug("This will make a simulation box for GCMC/JAWS-1")
  
   # Create a box around the solute and pad it with two Angstromgs
  box = pdb.getBox()
  box["origin"] = box["origin"] - padding
  box["len"] = box["len"] + 2.0*padding
  
  # Save it to disc
  simulationobjects.write_box(filename,box)

if __name__ == "__main__":

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to make a PDB-file with box coordinates covering a solute molecules")
  parser.add_argument('-s','--solute',help="the name of the PDB-file containing the solute.")
  parser.add_argument('-p','--padding',type=float,help="the padding in A,default=2",default=2.0)
  parser.add_argument('-o','--out',help="the name of the box PDB-file",default="gcmc_box.pdb")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("make_gcmcbox_py.log")

  if args.solute is None : 
    print "Nothing to do! Exiting."
  
  make_gcmcbox(args.solute,args.out,args.padding)
