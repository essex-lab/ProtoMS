# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to split a multi pdb file into individual files

This program divides a multi pdb file in individual files.
The individual pdbs must be separated by a line starting with END.

Initially thought to divide the multi-snapshot pdb file produced as ProtoMS output
into individual snapshot pdb files.
Required for an easy visualization of GCMC and JAWS stage 1 results.

Can be executed from the command line as a stand-alone program
"""

import os
import simulationobjects

def _get_prefix(filename) :
  """
  Remove extension (including period from a filename)

  Parameters
  ----------
  filename : string
    the filename to modify

  Returns
  -------
  string
    the filename without extension
  """
  h,t = os.path.splitext(filename)
  return h



def divide_print(common,prefix) :

  """
  Divide a pdb file in individual pdbfiles.
  Division is controled by lines starting with END

  Parameters
  ----------
  common : string
    the name of the initial pdb file
  prefix : string
    the prefix for output pdb files

  Returns
  -------
  """

  pdbset = simulationobjects.PDBSet()
  pdbset.read(filename=common)
  filenames = [prefix+"%d.pdb"%(i+1) for i,eachpdb in enumerate(pdbset.pdbs)]
  pdbset.write(filenames=filenames)

  return

def get_arg_parser ():
  import argparse

  
  # Setup a parser of the command-line arguments 
  parser = argparse.ArgumentParser(description="Split your multi pdb file into individual files")
  parser.add_argument('-i','--input',help="The name of your multi pdb file. Default = all.pdb",default="all.pdb")
  parser.add_argument('-o','--output',help="The basename of your individual pdb files. Default = snapshot_",default="snapshot_")
  parser.add_argument('-p','--path',help="Where the input should be found and the output printed. Default = ./",default="./")
  return parser


if __name__ == '__main__' :

  args = get_arg_parser().parse_args()

  # Make sure that the input file exists and get the prefix
  commonprefix = _get_prefix(args.input)

  commonfile = os.path.join(args.path,commonprefix+".pdb")

  if os.path.isfile(commonfile) :
    "The file %s will be divided in individual pdbs"%commonfile
  else :
    msg = "pdb file %s.pdb could not be found"%commonprefix
    print msg
    raise simulationobjects.SetupError(msg)

  # Get the prefix of the output files
  outprefix = _get_prefix(args.output)

  outpathprefix = os.path.join(args.path,outprefix)

  # Divide input into individual pdb objects
  divide_print(commonfile,outpathprefix)








