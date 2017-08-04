#!/usr/bin/env python2.7
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
  print 'Making GCMC box' 
   # Create a box around the solute and pad it with two Angstromgs
  box = pdb.getBox()
  box["origin"] = box["origin"] - padding
  box["len"] = box["len"] + 2.0*padding

  volume = volume_box(box["len"])
  print_bequil(volume)


  # Save it to disc
  simulationobjects.write_box(filename,box)


def make_gcmcsphere(pdb,filename,padding=2.0) :
  """
  Make a GCMC/JAWS-1 simulation sphere around a PDB-structure
  
  Parameters
  ----------
  pdb : PDBFile object
    the PDB structure
  filename : string
    the name of the output file
  padding : float, optional
    the amount of extra space around the ligand to add
  """
  
  logger.debug("Running make_gcmsphere with arguments: ")
  logger.debug("\tpdb      = %s"%pdb) 
  logger.debug("\tfilename = %s"%filename) 
  logger.debug("\tpadding  = %d"%padding) 
  logger.debug("This will make a simulation box for GCMC/JAWS-1")
  print 'Making GCMC sphere' 

 # finds the center of the GCMC sphere
  center,radius = pdb.getSphere()
  radius = radius + padding
  spherevol = volume_sphere(radius)
  print_bequil(spherevol)

  simulationobjects.write_sphere(filename,center,radius)

# print an atom with the centroid
#  simulationobjects.write_box(filename,box)

def volume_box(boxlen):
  return boxlen[0]*boxlen[1]*boxlen[2]

def calc_radius(pdb, center):
  print pdb
  print pdb.residue[0].coords

def volume_sphere(radius):
  return (4./3.)*np.pi*(radius**3.)

def print_bequil(volume):
  betamu = -10.47
  bequil = betamu+np.log(volume/30.0)
  print "Volume of GCMC region:", np.round(volume,2)
  print "Bequil:", np.round(bequil,2)

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to make a PDB-file with box coordinates covering a solute molecules")
  parser.add_argument('-s','--solute',help="the name of the PDB-file containing the solute.")
  parser.add_argument('-p','--padding',type=float,help="the padding of box or radius of sphere in A,default=2",default=2.0)
  parser.add_argument('-o','--out',help="the name of the box PDB-file",default="gcmc_box.pdb")
  parser.add_argument('--sphere',help="If flag given, a gcmc sphere will be defined",action='store_true')
  parser.add_argument('-b','--box',nargs='+',help="Either the centre of the box (x,y,z), or the centre of box AND length (x,y,z,x,y,z). If the centre is specified and the length isn't, twice the 'padding' will be the lengths of a cubic box.",default=None)
  return parser
  

if __name__ == "__main__":

  args = get_arg_parser().parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("make_gcmcbox_py.log")

  if args.solute is None and args.box is None : 
    print "Nothing to do! Exiting."  

  if args.sphere == True:
      pdbobj = simulationobjects.PDBFile()
      pdbobj.read(args.solute)
      make_gcmcsphere(pdbobj,'gcmc_sphere.pdb',args.padding)    
  else:
    if args.box is None :
      pdbobj = simulationobjects.PDBFile()
      pdbobj.read(args.solute)
      make_gcmcbox(pdbobj,args.out,args.padding)
    elif len(args.box) == 3:
      box = {"center":np.array([float(args.box[0]),float(args.box[1]),float(args.box[2])]),"len":np.array([args.padding*2]*3)}
      volume = volume_box(box["len"])
      print_bequil(volume)
      simulationobjects.write_box(args.out,box)
    elif len(args.box) == 6:
      box = {"center":np.array([float(args.box[0]),float(args.box[1]),float(args.box[2])]),"len":np.array([float(args.box[3]),float(args.box[4]),float(args.box[5])])}
      volume = volume_box(box["len"])
      print_bequil(volume)
      simulationobjects.write_box(args.out,box)
    else : 
      print "\nError with 'box' arguement. Please specify either three arguements for the centre of the box, or six arguements for the centre of the box AND the lengths of the sides.\n"
  
