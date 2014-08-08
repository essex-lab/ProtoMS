# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routines to convert water molecules to water models from a PDB file.

This module defines a single public function
convertwater
    
Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np
import simulationobjects 

logger = logging.getLogger('protoms')

def _translatetemplate(solvents,watresname,wattemplate,watatomnames):
  """ 
  Translates an ideal water model geometry, such as tip4p, to match the location of the water molecules in a residue object. Original hydrogen positions are over-written.

  Parameters
  ----------        
  solvents : dictionary of Residue objects 
  	  the residue object that contains the water oxygen positions
  watresname : string 
      the name of the water model that will be used
  wattemplate : numpy array
      the 3D location of the atoms (i.e. the geometry) of the ideal water model
  watatomnames : list of strings 
      the atom names of the water model as they will appear in the PDB file

  Returns
  -------
  dictionary of Residue objects
      a set of water molecules that have the 3D structure of an ideal water geometry
  """
  new_solvents = {}
  for sol in solvents:
      oxy_coord = solvents[sol].atoms[0].coords			# Assuming the water oxygen is the first entry in the solvent residue.
      translate = oxy_coord - wattemplate[0]					
      wat_new = wattemplate + translate													# Shifting the model water to match the oxygens position.
      new_solvents[sol] = simulationobjects.Residue(name=watresname,index=sol)			# Creating a new pdb object that contains the water.
      for ind in range(len(watatomnames)):
          newatom = simulationobjects.Atom(index=ind,name=watatomnames[ind],resindex=sol,resname=watresname,coords=wat_new[ind])
          new_solvents[sol].addAtom(atom=newatom)
  return (new_solvents)


def convertwater(pdb_in,watermodel):
  """ 
  Converts water in a pdb object to ideal water geometries of an input model
  
  The protein object is not modified by this routine, but the Residue and Atom objects are.  
    
  Parameters
  ----------        
  pdb_in : PDBFile 
      the PDB object containing the water molecules that will be converted
  watermodel : string 
      the name of the water model that will be used in the transformtation, e.g. tip4p, t3p.
  
  Returns
  -------
  PDBFile
      a pdb file whose solvent elements have the geometry of the desired solvent model
  """
  
  logger.debug("Running convertwater with arguments: ")
  logger.debug("\tpdb_in     = %s"%pdb_in) 
  logger.debug("\twatermodel = %s"%watermodel) 
  logger.debug("This will change the water molecule in the pdb file to match the water model")
  
  pdb_out = pdb_in.copy()
  solvents = pdb_in.solvents
  # Ideal water geometries:
  t4p_model = np.array([[-6.444, -5.581, -1.154],[-7.257, -5.869, -0.739],[-6.515, -4.627, -1.183],[-6.557, -5.495, -1.105]]) 
  t4p_names = ["O00","H01","H02","M03"]
  t3p_model = np.array([[7.011, 5.276, 13.906],[6.313, 5.646, 13.365],[6.666,  4.432, 14.197]]) 
  t3p_names = ["O00","H01","H02"]
  if watermodel.upper() in ["T4P","TIP4P","TP4"]: 
      pdb_out.solvents = _translatetemplate(solvents,"T4P",t4p_model,t4p_names)
  elif watermodel.upper() in ["T3P","TIP3P","TP3"]:
      pdb_out.solvents = _translatetemplate(solvents,"T3P",t3p_model,t3p_names)
  else:
      print "Error in convertwater.py: water model name not recognised. Please check spelling matches known list or add new water model to function." 
  return pdb_out

if __name__ == "__main__":

    import argparse
    
    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to convert water molecules - with or without hydrogens - in a pdb file to simulation models, such as tip4p. Currently ignores original hydrogen positions.")
    parser.add_argument('-p','--pdb',help="the PDF-file containing the waters to be transformed")
    parser.add_argument('-o','--out',help="the output PDB-file",default="convertedwater.pdb")
    parser.add_argument('-m','--model',help="the water model,default=tip4p",default="tip4p")
    args = parser.parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("convertwater_py.log")

    pdb_in = simulationobjects.PDBFile(filename=args.pdb)
    pdb_out = convertwater(pdb_in,args.model)
    pdb_out.write(args.out)

