# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to generate a set of n molecules contained within given box dimenstions
or to redistribute molecules given in a pdb object or file within the given box dimensions

This program takes the box dimensions as a list of strings convertible to floats,
in the following order :
  origin in x | origin in y | origin in z | length in x | length in y | length in z
And either the number of molecules as a string convertible to integer,
or the pdb, which can be provided as a string with the name of the file, or an object.

It produces a pdb file with the randomly distributed molecules with their oxygen atoms
within the box limits.

Initially thought to generate JAWS or GCMC molecules

Can be executed from the command line as a stand-alone program
"""

import simulationobjects
import logging
import convertwater
import numpy as np

logger = logging.getLogger('protoms')

def create_res(atnames=["O00"],resname="sol",positions=[np.array([0.0,0.0,0.0])],resind=1) :
  """
  Create a residue object

  Parameters
  ----------
  atnames : list of strings
    a list containing the names of all
    the atoms in the residue
  resname : string
    the name of the residue
  positions : list of numpy arrays
    a list with the positions of all
    the atoms in the residue
  
  Returns
  -------
  Residue object
    the residue built with the given atoms
  """
  resobj = simulationobjects.Residue()
  resobj.name = resname
  resobj.index = resind
  for ind,atom in enumerate(atnames) :
    atobj = simulationobjects.Atom()
    atobj.name = atom
    atobj.resname = resname
    atobj.index = ind
    atobj.resindex = resind
    try :
      atobj.coords = positions[ind]
    except :
      atobj.coords = np.array([0.0,0.0,0.0])
    resobj.addAtom(atobj)
  return resobj
  



def distribute_particles(box, particles, watermodel="t4p",resname="WA1",partnumb=None) :
  """
  Randomly distribute molecules in a box
  
  Parameters
  ----------
  box : a list or dictionary
    the origin and length of the box
    where the molecules will be placed
  particles : string or PDB object
    the number of waters to include in the box
    or the pdb object or filename with the molecules
    to include in the box
  watermodel : string, optional
    (only used when particles is a number)
    for instance, "t4p", "tip4p", "t3p" "tip3p",
    the water model for the generated waters
  partnumb : string,optional
    (only used when particles is a file)
    the number of particles. If not specified it is set
    to be as many as there are in the file

  Returns
  -------
  pdb object
    the pdb object with the molecules
    distributed randomly in the box
  """

  logger.debug("Running distribute_particles with arguments: ")
  if isinstance(box,dict) :
    tmpstr = "\tbox = "
    for key in box:
      tmpstr = tmpstr + "%s %s, "%(key,str(box[key]))
    logger.debug(tmpstr.strip(" ,"))
  else :
    logger.debug("\tbox = %s"%" ".join(box)) 
  if isinstance(particles,str) :
    logger.debug("\tparticles = %s"%particles)
  else :
    logger.debug("\tparticles %s"%particles.name)
  logger.debug("\twatermodel %s"%watermodel)
  logger.debug("\tpartnumb %s"%partnumb)
  logger.debug("This will distribute the molecules given in 'particles' randomly in the box given in box. If a number is given in 'particles', it assumes them to be waters.")

  if isinstance(box,list) :
    try :
      box = [float(val) for val in box]
    except :
      raise simulationobjects.SetupError("The box dimensions %s could not be correctly interpreted"%box)
    orig = box[:3]
    length = box[3:]
  elif isinstance(box,dict) :
    orig = box['origin']
    length = box['len']

  if isinstance(particles,str) and particles.isdigit():
    watnumb = int(particles)
    particles = simulationobjects.PDBFile()
    for i in range(1,watnumb+1) :
      pad = 0.3                 # Water oxygens are at least this far away (in Angs.) from the edge of the box
      oxpos = np.random.rand(len(length))*np.array(length-2*pad) + np.array(orig+pad)
      particles.solvents[i] = create_res(resname=resname,positions=[oxpos])
    particles = convertwater.convertwater(particles,watermodel,"y",watresname=resname)

  else :
    if type(particles) is str :
      try :
        pdbobj = simulationobjects.PDBFile()
        pdbobj.read(particles)
      except :
        raise simulationobjects.SetupError("The pdb file %s could not be found"%particles)
      particles = pdbobj
    if particles.residues :
      parts = particles.residues
    elif particles.solvents : 
      parts = particles.solvents
    else :
      raise simulationobjects.SetupError("No molecule could be found in %s"% particles)
    if partnumb is not None :
      particles = simulationobjects.PDBFile()
      for i in range(1,int(partnumb)+1) :
        allcoords = [myat.coords for myat in parts[parts.keys()[0]].atoms]
        allnames = [myat.name for myat in parts[parts.keys()[0]].atoms]
        particles.residues[i] = create_res(resname=parts[parts.keys()[0]].name,positions=allcoords,atnames=allnames,resind=i)
      parts = particles.residues    
    for ind,keypart in enumerate(parts) :
      parts[keypart].getCenter()
      displace = [coord_len*np.random.uniform() + orig[jnd] - parts[keypart].center[jnd] for jnd,coord_len in enumerate(length)]
      rotated_coords = convertwater.rotatesolute(np.array([myat.coords for myat in parts[keypart].atoms]),np.random.uniform(0,2*np.pi),np.random.uniform(0,2*np.pi),np.random.uniform(0,2*np.pi))
      for jnd,atom in enumerate(parts[keypart].atoms) :
        atom.coords = np.array([rotated_coords.item((jnd,i)) + displace[i] for i in range(3)])
      

  for ind,coord in enumerate(orig) :
    h_parts = particles.header.strip().split()
    particles.header = "%s %.4f %s %.4f " %(" ".join(h_parts[:ind]),coord," ".join(h_parts[ind:ind*2]),coord+length[ind])
  particles.header = "HEADER box %s\n"%particles.header
  return particles
    

if __name__ == "__main__":

  import argparse
  
  parser = argparse.ArgumentParser(description="Randomly distribute n molecules within box dimensions")
  parser.add_argument('-b','--box',nargs=6,help="Dimensions of the box. Six arguments expected: origin (x,y,z) & length (x,y,z)")
  parser.add_argument('-m','--molecules',help="Molecules to distribute in the box. Either the number of waters or a pdb file containing all of them")
  parser.add_argument('-o','--outfile',help="Name of the pdb file to write the molecules to. Default='ghostmolecules.pdb'",default='ghostmolecules.pdb')
  parser.add_argument('--model', help="Water model. Used when only the amount of waters is specified. Options: 't4p','t3p'. Default='t4p'",default='t4p')
  parser.add_argument('--resname',help="Residue name of the molecules writen to output. Default='WAT'",default='WAT')
  parser.add_argument('--number',help="Required number of molecules when it differs from the number of residues in the file.",default=None)
  parser.add_argument('--setupseed',help="Optional random number seed for generation of water coordinates.",default=None,type=int)
  args = parser.parse_args()

  logger = simulationobjects.setup_logger("distribute_waters_py.log")
  if args.setupseed is not None:
    logger.debug("Setup seed = %d"%args.setupseed)
  np.random.seed(args.setupseed)

  outobj = distribute_particles(args.box,args.molecules,args.model,args.resname,args.number)
  outobj.write(filename=args.outfile)
  print "\nMolecules printed in %s"%args.outfile

