# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routines to estimate the free volume within certain box.

It also estimates an upperlimit on how many solutes are needed to
fill the cavity volume.

This set of routines are grid based, where distance between
grid points and the approximation to the "Van de Waals radius"
of my protein atoms can be modified.
    
Can be executed from the command line as a stand-alone program
"""

import logging
import numpy as np
import simulationobjects
from distribute_waters import create_res
logger = logging.getLogger('protoms')



def load_pdb(pdbfile) :
  """
  Function which attempts at loading a pdbfile when
  the argument given is a string

  Parameters
  ----------
  pdbfile: string or PDBFile object
    the string is attempted to be opened as
    a PDBFile object. If the argument passed
    is not a string, it is returned without changes

  Returns
  -------
  PDBFile object
    the object corresponding to the filename read
    of the argument passed as input wihout any changes
  """

  if isinstance(pdbfile,str) :
    try :
      pdbfile = simulationobjects.PDBFile(filename=pdbfile)
    except :
      raise simulationobjects.SetupError("The pdbfile %s could not be oppened."%pdbfile)

  return pdbfile
    


def set_grid(box,grid_dist=1.0) :
  """
  Function that generates the coordinates of
  the grid points

  Parameters
  ----------
  box: dict
    a dictinary or numpy arrays
    with at least two entries
    with dimensions for box length and
    box origin and/or center
  grid_dist: float, optional
    the distance between two grid points
    in angstroms

  Returns
  -------
  numpy array
    the xyz coordinates of the grid points
  """

  if not isinstance(box,dict) :
    raise simulationobjects.SetupError("The box information passed has not the expected format. Dict required.")
  elif (not 'origin' in box.keys() and not 'center' in box.keys()) or (not 'len' in box.keys()) :
    raise simulationobjects.SetupError("Missing information on the box where the grid should be applied.")

  for key in box :
    if not isinstance(box[key], np.ndarray) :
      raise simulationobjects.SetupError("Numpy array expected for box %s"%key)

  if not grid_dist > 0 :
    raise simulationobject.SetupError("The grid distance must be positive")

  if not 'origin' in box.keys() :
    box['origin'] = box['center'] - 0.5*box['len']
 

  # Generating list of lists with all coordinates in each of the axes (i.e. all x coordinates in coords[0])

  coords = [np.arange(orig_coord+grid_dist/2,box['len'][ind]+orig_coord-grid_dist/2+grid_dist/100,grid_dist) for ind,orig_coord in enumerate(box['origin'])]

  grid_points = np.array([[i,j,k] for i in coords[0] for j in coords[1] for k in coords[2]])

  return grid_points



def clear_grid(grid,at_vdw=1.5,pdbfiles=None) :
  """
  Function which calculates the set of grid points
  which do not overlap with any atom of the system

  Parameters
  ----------
  grid: numpy array
    a numpy array with all points of the grid which
    fit withint the box dimensions
  at_vdw: float, optional
    the estimation of the VdW radius of the atoms in
    the pdb files. Since all calculations are approximate
    the same VdW radius is used for every atom
  pdbfiles: list, optional
    a list of PDBFiles containing any ligand or protein
    which space should be excluded from the grid

  Return
  ------
  numpy array
    the cleared grid
  list of numpy arrays
    list of arrays with the grid occupied by each of
    the pdb files included in pdbfiles
  """

  cleared_grid = grid.tolist()

  if pdbfiles :
    pdb_grids = [[] for i,pdb in enumerate(pdbfiles)]
    for ind,pdbfile in enumerate(pdbfiles) :
      pdbfile = load_pdb(pdbfile)
      for res in pdbfile.residues :
        for atom in pdbfile.residues[res].atoms :
           at_dist = at_vdw
           atom.getElement()
           if atom.element.lower() in 'h' or atom.element.lower() in 'm' : at_dist = at_vdw*0.5
           for point in cleared_grid :
             dist = (np.sum([(at_coord-point[jnd])**2 for jnd,at_coord in enumerate(atom.coords)]))**0.5
             if dist < at_dist :
               pdb_grids[ind].append(point)
           cleared_grid = [point for point in cleared_grid if not point in pdb_grids[ind]]
    pdb_grids[ind] = np.array(pdb_grids[ind])

  return np.array(cleared_grid),pdb_grids



def grid_to_pdb(grid,pdbname="cavity_grid.pdb") :
  """
  Function to generate a pdbfile out of 
  the grid

  Parameters
  ----------
  grid: numpy array
    numpy array containing the grid points
  pdbname: string,optional
    name of the pdb file where the grid
    points will be printed

  Returns
  -------
  None
    It prints the pdbfile in pdbname
  """

  pdbobj = simulationobjects.PDBFile()
  for ind,point in enumerate(grid) :
    pdbobj.solvents[ind+1] = create_res(positions=[point],resind=ind+1,resname="GRD")

  if len(grid) > 9999 :
    print "\nOwing to huge number of grid 'atoms', %s does not fulfil pdb standards. It's visualization will be problematic."%pdbname

  pdbobj.write(filename=pdbname)



def cavity_volume(box,grid_dist,vdw_rad=1.5,cav_pdbs=None,sol_pdbs=None) :
  """
  Function which generates a grid on the cavity specified,
  estimates a corresponding free space volume of the cavity,
  as well as the number of solute copies to fill the cavity.
  
  Parameters
  ----------
  box: dict
    the dictionary containing the origin and length of the
    box comprising the cavity
  grid_dist: float
    the length of one side of the cubic boxel
    which volume each grid point represents. In
    angstroms
  vdw_rad: float, optional
    the estimation of the VdW radius of the atoms in
    the pdb files. Since all calculations are approximate
    the same VdW radius is used for every atom
  cav_pdbs: list, optional
    list of filenames or PDBFiles of the molecules meant
    to define the cavity
  sol_pdbs: list, optional
    list of filenames or PDBFiles of the solutes for which
    the number of copies needed to fill the cavity will be
    calculated

  Returns
  -------
  float
    the volume in angstroms, of the cavity, estimated with
    the grid
  integer
    the number of copies of each solute needed to fill the cavity
    if all solutes are to be in the same proportion
  list of integers
    the number of copies required to fill the cavity with each of
    the solutes individually
  """

  logger.debug("Running cavity_volume with arguments: ")
  logger.debug("\tbox      = %s"%box) 
  logger.debug("\tgrid_dist = %.2f"%grid_dist) 
  logger.debug("\tvdw_rad  = %.2f"%vdw_rad) 
  logger.debug("\tcav_pdbs = %s"%cav_pdbs) 
  logger.debug("\tsol_pdbs  = %s"%sol_pdbs)
  logger.debug("This will make a simulation box for GCMC/JAWS-1")

  grid = set_grid(box,grid_dist)

  if cav_pdbs : grid,cav_grids = clear_grid(grid,pdbfiles=cav_pdbs)

  pdbname="cavity_grid.pdb"
  grid_to_pdb(grid,pdbname=pdbname)

  logger.info("Grid pdb for %s printed in %s"%("cavity",pdbname))

  cavity_vol = len(grid)*grid_dist**3

  if sol_pdbs :
    grid_vols = []
    for ind,sol_pdb in enumerate(sol_pdbs) :
      sol_obj = load_pdb(sol_pdb)
      box_sol = sol_obj.getBox()
      box_sol["origin"] = box_sol["origin"] - vdw_rad
      box_sol["len"] = box_sol["len"] + 2.0*vdw_rad
      solbox_grid = set_grid(box_sol,grid_dist)
      null_grid,sol_grid = clear_grid(solbox_grid,pdbfiles=[sol_obj])
      gridpdb_name = "pdb%d_grid.pdb"%(ind+1)
      grid_to_pdb(sol_grid[0],pdbname=gridpdb_name)
      grid_vols.append(len(sol_grid[0]))

      logger.info("Grid pdb for %s printed in %s"%("solute %d"%(ind+1),gridpdb_name))

    each_copies = [round(len(grid)/grid_vol) for grid_vol in grid_vols]
    all_copies = round(len(grid)/sum(grid_vols)) 

    return cavity_vol,all_copies,each_copies

  else :

    return cavity_vol,None,[None]



if __name__ == '__main__' :
 
  import argparse

  parser = argparse.ArgumentParser(description="Program to estimate the volume of a cavity and generate an upper limit estimate of the number of copies of solutes required to fill the cavity")
  parser.add_argument('-o','--boxorigin',nargs='+',help="The coordinates of the bottom left down corner of the box")
  parser.add_argument('-l','--boxlength',nargs='+',help="The coordinates of the length of the box")
  parser.add_argument('-d','--distance',default=1,help="The distance between two grid points.Default=1(A)")
  parser.add_argument('-w','--vdwrad',default=1.5,help="The estimated Van de Waals radius for most atoms. It is halved for hydrogen and the lone pair in T4P. Default=1.5(A)")
  parser.add_argument('-c','--cavpdbs',nargs='+',default=None,help="The pdb files of the molecules definining the cavity")
  parser.add_argument('-s','--solpdbs',nargs='+',default=None,help="The pdb files of the solutes that we wish to fit in the cavity")
  parser.add_argument('-v','--volume',action="store_true",default=False,help="Whether to print the volume of the cavity. Default=False")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("cavity_volume.py.log")

  box_dict={'origin': np.array(args.boxorigin,dtype=float),'len': np.array(args.boxlength,dtype=float)}

  try :
    vdw_rad = float(args.vdwrad)
  except :
    raise simulationobjects.SetupError("VdW radius specified could not be read as float")

  try :
    grid_dist = float(args.distance)
  except :
    raise simulationobjects.SetupError("Distance specified could not be read as float")

  cav_vol,all_copies,each_copies = cavity_volume(box_dict,grid_dist,vdw_rad,args.cavpdbs,args.solpdbs)

  # Printing results below

  if args.volume : logger.info("Estimated cavity volume: %.2f A^3" %cav_vol)

  if args.solpdbs :

    logger.info("\nUpper-estimate: %d copies of each solute to fill cavity volume" %all_copies)

    for ind,solute in enumerate(args.solpdbs) :
      logger.info("\nUpper-estimate: %d copies of solute %s to fill cavity volume" %(each_copies[ind],solute))






