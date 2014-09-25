# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to generate a set of n waters contained within given box dimenstions
or to redistribute waters given in a pdb object or file within the given box dimensions

This program takes the box dimensions as a list of strings convertible to floats,
in the following order :
  origin in x | origin in y | origin in z | length in x | length in y | length in z
And either the number of waters as a string convertible to integer,
or the pdb, which can be provided as a string with the name of the file, or an object.

It produces a pdb file with the randomly distributed waters with their oxigen atoms
within the box limits.

Initially thought to generate JAWS or GCMC waters

Can be executed from the command line as a stand-alone program
"""

import simulationobjects
import random
import convertwater
import numpy as np



def create_res(atnames=["O00"],resname="sol",positions=[np.array([0.0,0.0,0.0])]) :
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
  for ind,atom in enumerate(atnames) :
    atobj = simulationobjects.Atom()
    atobj.name = atom
    try :
      atobj.coords = positions[ind]
    except :
      atobj.coords = np.array([0.0,0.0,0.0])
    resobj.addAtom(atobj)
  return resobj
  



def distribute_waters(box, waters, watermodel="t4p", out="ghostwaters.pdb",resname="WAT") :
  """
  Randomly distribute waters in a box
  
  Parameters
  ----------
  box : a list or dictionary
    the origin and lenght of the box
    where the waters will be placed
  waters : string or PDB object
    the number of waters to include in the box
    or the pdb object or filename with the waters
    to include in the box
  watermodel : string, optional
    (only used when waters is a string)
    either "t4p" or "t3p"
    the water model for the generated waters
  out : string, optional
    the name of the pdb file to write the waters to

  Returns
  -------
  string
    the pdb file name of the file with the waters
    distributed randomly in the box
  """

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

  if isinstance(waters,str) and waters.isdigit():
    watnumb = int(waters)
    waters = simulationobjects.PDBFile()
    for i in range(1,watnumb+1) :
      oxpos = np.array([coord_len*random.random()+orig[jnd] for jnd,coord_len in enumerate(length)])
      waters.solvents[i] = create_res(resname=resname,positions=[oxpos])
    waters = convertwater.convertwater(waters,watermodel,"y",watresname=resname)

  else :
    if type(waters) is str :
      try :
        pdbobj = simulationobjects.PDBFile()
        pdbobj.read(waters)
      except :
        raise simulationobjects.SetupError("The pdb file %s could not be found"%waters)
      waters = pdbobj
    allwats = len(waters.solvents)
    for ind,solv in enumerate(waters.solvents) :
      displace = [coord_len*random.random() + orig[jnd] - waters.solvents[solv].atoms[0].coords[jnd] for jnd,coord_len in enumerate(length)]
      for atom in waters.solvents[solv].atoms :
        atom.coords = np.array([atom.coords[i] + displace[i] for i in range(3)])

  for ind,coord in enumerate(orig) :
    h_parts = waters.header.strip().split()
    waters.header = "%s %.4f %s %.4f " %(" ".join(h_parts[:ind]),coord," ".join(h_parts[ind:ind*2]),coord+length[ind])
  waters.header = "HEADER box %s\n"%waters.header
  waters.write(filename=out)
  return out
    

if __name__ == "__main__":

  import argparse
  
  parser = argparse.ArgumentParser(description="Randomly distribute n waters within box dimensions")
  parser.add_argument('-b','--box',nargs='+',help="Dimensions of the box. Six arguments expected: origin (x,y,z) & length (x,y,z)")
  parser.add_argument('-w','--waters',help="Waters to distribute in the box. Either the number of waters or a pdb file containing all of them")
  parser.add_argument('-o','--outfile',help="Name of the pdb file to write the waters to. Default='ghostwaters.pdb'",default='ghostwaters.pdb')
  parser.add_argument('--model', help="Water model. Used when only the amount of waters is specified. Options: 't4p','t3p'. Default='t4p'",default='t4p')
  parser.add_argument('--resname',help="Residue name of the waters writen to output. Default='WAT'",default='WAT')
  args = parser.parse_args()

  if len(args.box) < 6 :
    raise simulationobjects.SetupError("Not enough information regarding box dimensions: %s"%args.box)

  outfile = distribute_waters(args.box,args.waters,args.model,args.outfile,args.resname)
  print "\nWaters printed in %s"%outfile

