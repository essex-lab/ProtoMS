# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routine to split JAWS-1 water into PDBs for JAWS-2

This module defines a single public function:
split_water

Can be executed from the command line as a stand-alone program
"""

import logging
import copy

import numpy as np

import simulationobjects

logger = logging.getLogger('protoms')

def split_waters(waters) :
  """
  Split waters in a PDB file to many PDB-files for JAWS-2
  
  Parameters
  ----------
  waters : string or PDBFile object
    the PDB structure

  Returns
  -------
  PDBSet
    single waters
  PDBSet
    excluding single waters
  """
  
  logger.debug("Running split_waters with arguments: ")
  logger.debug("\twaters = %s"%waters) 
  logger.debug("This will make PDB files suitable for JAWS-2")
  
  if isinstance(waters,basestring) :
    waters = simulationobjects.PDBFile(filename=waters)

  single_waters = simulationobjects.PDBSet()
  single_waters.from_residues(waters)

  other_waters = simulationobjects.PDBSet()
  for soli,sol in waters.solvents.iteritems() :
    others = copy.deepcopy(waters)
    del others.solvents[soli]
    other_waters.pdbs.append(others)
    others.header = ""
    # New to change the residue name of the solvent molecules, this is a hack for now
    for soli2,sol2 in others.solvents.iteritems() :
      if len(sol2.atoms) == 3 :
        resnam = "t3p"
      elif len(sol2.atoms) == 4 :
        resnam = "t4p"
      for atom in sol2.atoms : atom.resname = resnam

  return single_waters,other_waters

if __name__ == "__main__":

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to split JAWS-1 waters to a number of PDB-files for JAWS-2")
  parser.add_argument('-w','--waters',help="the name of the PDB-file containing the waters.")
  parser.add_argument('-o','--out',help="the prefix of the output PDB-files",default="")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("split_jawswater_py.log")

  if args.waters is None : 
    print "Nothing to do! Exiting."
  
  single_waters,other_waters = split_waters(args.waters)

  single_waters.write([args.out+"wat%d.pdb"%(i+1) for i in range(len(single_waters.pdbs))])
  other_waters.write([args.out+"not%d.pdb"%(i+1) for i in range(len(single_waters.pdbs))])
