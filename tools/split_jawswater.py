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
  for soli,sol in waters.residues.iteritems() :
    others = copy.deepcopy(waters)
    del others.residues[soli]
    other_waters.pdbs.append(others)
    others.header = ""
    # New to change the residue name of the molecules, this is a hack for now
    for soli2,sol2 in others.residues.iteritems() :
      if len(sol2.atoms) == 3 :
        resnam = "t3p"
      elif len(sol2.atoms) == 4 :
        resnam = "t4p"
      for atom in sol2.atoms : atom.resname = resnam

  return single_waters,other_waters

def set_jaws2_box(water_files,length=3) :
  """
  Sets the header of the pdbs containing one single
  molecule to define a box around the first atom

  Parameters:
  ----------
  water_files : PDBSet
    a pdb set containing the pdb files with one molecule
  length : int, optional
    the dimensions of the box around the molecule

  Returns:
  --------
  PDBSet
    the same input PDBSet, with changed headers
  """
  for count,watobj in enumerate(water_files.pdbs) :
    boxcoords= np.zeros(3)
    for i,coord in enumerate(watobj.residues[watobj.residues.keys()[0]].atoms[0].coords) :
      boxcoords[i] = coord-1.5
      boxcoords = np.append(boxcoords,coord+1.5)
    watobj.header = watobj.header + "REMARK box"
    for coord in boxcoords :
      watobj.header = watobj.header + " %.3f"%coord
    watobj.header = watobj.header + "\n"
  return water_files

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to split JAWS-1 waters to a number of PDB-files for JAWS-2")
  parser.add_argument('-w','--waters',help="the name of the PDB-file containing the waters.")
  parser.add_argument('-o','--out',help="the prefix of the output PDB-files",default="")
  parser.add_argument('--jaws2box',action='store_true',help="whether to apply a header box for jaws2 to the pdb files of individual waters",default=False)
  return parser

if __name__ == "__main__":

  args = get_arg_parser().parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("split_jawswater_py.log")

  if args.waters is None : 
    print "Nothing to do! Exiting."
  
  single_waters,other_waters = split_waters(args.waters)

  if args.jaws2box :
    set_jaws2_box(single_waters)

  single_waters.write([args.out+"wat%d.pdb"%(i+1) for i in range(len(single_waters.pdbs))])
  other_waters.write([args.out+"not%d.pdb"%(i+1) for i in range(len(single_waters.pdbs))])










