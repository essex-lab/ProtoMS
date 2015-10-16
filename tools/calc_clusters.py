# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routines to cluster molecules from simulations.

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance

import simulationobjects

logger = logging.getLogger('protoms')


def get_coords(pdbfiles,molname,atomname):
  """
  Finds the coordinates of a specified molecule name or a named atom from supplied PDB files.
  
  Parameters
  ----------
  pdbfiles : PDBSet
    the PDB files
  residue : string
    the residue to extract
  atom : string
    the name of the atom to extract (optional)

  Returns
  -------
  NumpyArray of coordinates
    the x,y,z coordinates of the specified molecule or atom
  Integer
     the number of molecules or atoms matching the names.
  """
  # Extract coordinates from PDB-files
  mol_xyz = []
  atm_xyz = []
  molsizes = [] 
  mol_found = 0
  atm_found = 0
  for pdb in pdbfiles.pdbs :
    for i,res in pdb.residues.iteritems() : # First, trying to locate the molecule in the list of residues (which is everything but water).
      if res.name.lower() != molname.lower() : continue
      mol_found = mol_found + 1
      for atom in res.atoms :
        if atomname is not None:
          if atom.name.strip().lower() == atomname.lower() :
            atm_xyz.append(atom.coords)
            atm_found = atm_found + 1
        else:
            mol_xyz.append(atom.coords)
            molsizes.append(len(res.atoms))
    for i,sol in pdb.solvents.iteritems() : # Second, trying to locate the molecule in the list of solvents.
      if sol.name.lower() != molname.lower() : continue
      mol_found = mol_found + 1
      for iatom in sol.atoms :
        if atomname is not None:
          if iatom.name.strip().lower() == atomname.lower() :
            atm_xyz.append(iatom.coords)
            atm_found = atm_found + 1
        else:
            mol_xyz.append(iatom.coords)
            molsizes.append(len(sol.atoms))
  if not mol_xyz :
    raise simulationobjects.SetupError("No molecule with residue name %s has been found."%molname)
  if atomname is not None :
    conformations = np.array(atm_xyz)
    return conformations, atm_found
  else :
    mols = np.array(molsizes)
    if np.std(mols) > 1E-6 : 
      print "\nInconistent number of atoms for molecule %s. They must all be the same molecule to be able to compute RMSDs. Aborting clustering.\n" % molname
      quit()
    conformations = np.array(mol_xyz).reshape(mol_found,molsizes[0],3)
    return conformations, mol_found

def cluster_coords(confs,clustmethod='average',cutoff=2.0):
  """
  Clusters the supplied conformations in a hierarchical fasion
  
  Parameters
  ----------
  conformations : Numpy array
    Array of the coordinates of the different conformations
  method : String
    the hierarchical clustering method
  cutoff : Float
    the RMSD cutoff that specifies whether conformations are in the same cluster.

  Returns
  -------
  NumpyArray 
    the coordinates of the cluster representatives
  NumpyArray
     the coordinates of the cluster centroids
  """
  # Calculate the pairwise RMSDs between the conformations.
  RMSDs = np.zeros((confs.shape[0], confs.shape[0]))
  for i in range(confs.shape[0]):
    x=confs-confs[i]
    y=np.sum(x**2,axis=2)
    RMSDs[i,:] = np.mean(y,axis=1)**(1./2)
  RMSDs = distance.squareform(RMSDs) 					# Correcting the formating for clustering.

  # The clustering bit.
  z = hierarchy.linkage(RMSDs,method=clustmethod,metric='euclidean')
  Clusts = hierarchy.fcluster(z,t=cutoff,criterion='distance')

  # Getting cluster representatives
  Occ = range(Clusts.max()) # Pre-assigning the occupancy (number in each cluster).
  ClustReps = np.zeros((Clusts.max(),confs.shape[1],3))	 # Pre-assigning the cluster representatives.
  ClustMeans = np.zeros((Clusts.max(),confs.shape[1],3)) # Pre-assigning the atomic position means of each cluster. These are NOT the same as the cluster representatives.
  MSD = np.zeros((confs.shape[1],Clusts.max())) # To hold the atomic mean deviations of the cluster representatives.
  for i in np.arange(1,Clusts.max()+1):	# Looping over the clusters and saving the conformations with the smallest distance from the atomic centroids.
    Occ[i-1] = np.sum(Clusts==i)
    if Occ[i-1]==1: # If there is only one molecule in the cluster, no need to average over a bunch as it is its own representative.
      ClustReps[i-1] = confs[np.where(Clusts==i)] # (Note that MSD is zero here, as in the pre-assignment).
      ClustMeans[i-1] = confs[np.where(Clusts==i)]
    elif Occ[i-1]==2: # If there are two in a cluster, the mean is equidistant from them both, so simply outputing the first in the cluster.
      InClusters = confs[np.where(Clusts==i)]
      ClustMeans[i-1] = InClusters.mean(axis=0)
      AtomicSrdDists = ((ClustMeans[i-1] - InClusters)**2).sum(axis=2) # The squared distance of each atom to centroid.
      MSD[:,i-1] = (AtomicSrdDists[0]).round(decimals=2) # = Ave(dist^2)
      ClustReps[i-1] = confs[np.where(Clusts==i)][0] # Choosing the first cluster member as the representative of the two. The Bfactor quantifies the uncertainty in this.
    else: # More than two in a cluster means that is it is meaningful to pick the closest to the mean as a representative.
      InClusters = confs[np.where(Clusts==i)]
      ClustMeans[i-1] = InClusters.mean(axis=0)  # The centroid of each atom of the cluster.
      AtomicSrdDists = ((ClustMeans[i-1] - InClusters)**2).sum(axis=2)	# The squared distance of each atom to centroid.
      SqrdDists = AtomicSrdDists.sum(axis=1) # The squared distance of each conformation to centroid.
      ClustReps[i-1] = InClusters[np.where(SqrdDists==SqrdDists.min())][0] # Picking the cluster representative closest to centroids. If there is a tie, pick the first.
      MSD[:,i-1] = (AtomicSrdDists[np.where(SqrdDists==SqrdDists.min())][0]).round(decimals=2)  # = Ave(dist^2).
  # Sorting in the order of cluster occupancy, from highest to lowest.
  #SortByOcc = sorted(zip(Occ,range(Clusts.max())),reverse=True)
  #SortByOcc = [y for (x,y) in SortByOcc]
  return ClustReps, np.array(Occ)

#
# Helper routines
#

def _GetMolTemplate(FileName, MolName):
	"""
        Extracts the pdb of the molecule of interest. This is then used as a template for printing out coordinates of that molecule.
        FileLines: a pdb file which is the output from a molecular simulation, and which may contain the coordinates of many other molecules.
        molname: the residue name of the molecule one wishes to extract theta values for, eg "wa1".
	"""
	FileLines = open(FileName,"r").readlines()
	molpdb = []
	for i in range(len(FileLines)): # Looping over all lines in pdb.
		line = FileLines[i]
		if line[0:4]=="ATOM" or line[0:6]=="HETATM":												
			if line[16:20].strip().lower() == MolName.lower(): # Finding the atoms that match the name of the molcule.
				molpdb.append(line)
		elif (line[0:3] == 'TER' or line[0:3] == 'END') and len(molpdb) > 0 :						
			molpdb.append('TER') # At the end of the first molecule, this program stops as 
			break # a complete pdb of the molecule has been extracted. 
	return molpdb				

def _printpdb(coords,molpdb,resid,theta,filename):
        """
        Prints the coordinates of a molecule in a pdb format, using a template pdb file. 
 	coords: the coordinates of the molcule one wishes to print in pdb format.
        molpdb: the template pdb file of the molecule.
        """
	for i in range(len(molpdb)):
		template = bytearray(molpdb[i]).replace('\n','')
		if i < (len(molpdb)-1):	# Altering only the coordinates, and leaving the TER line untouched.
			template[30:55]=' '*25	# Clearing the protoms junk.
			template[30:54]='{:>8}{:>8}{:>8}'.format(coords[i][0],  coords[i][1],  coords[i][2])
			template[22:26]='{:>4}'.format(resid) # Residue id by order of appearance.
			template[56:76]=' '*20	# Clearing the protoms junk.
			template[56:60]='{:<4}'.format(round(theta,2))	# Printing the molecules theta value in the occupancy column.
		print >> filename, template



if __name__ == "__main__":

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to cluster molecules")
  parser.add_argument('-f','--files',nargs="+",help="the input PDB-files")
  parser.add_argument('-o','--out',help="the name of the output pdb that contains the cluster representatives, default='clusters.pdb'",default="clusters.pdb")
  parser.add_argument('-m','--molecule',help="the name of the molecule to extract, default='wa1'",default="wa1")
  parser.add_argument('-a','--atom',help="the name of the atom to extract. If left blank, the entire molecule will be clustered, default=None",default=None)
  parser.add_argument('-t','--type',choices=["average","single","complete","weighted", "centroid"],help="the type of hierarchical clustering method, default='average'", default="average")
  parser.add_argument('-c','--cutoff',type=float,help="the distance cutoff that defines whether conformations belong to a cluster, default=2.0",default=2.0)
  parser.add_argument('--skip',type=int,help="the number of blocks to skip to calculate the clusters. default is 0. Skip must be greater or equal to 0",default=0)
  parser.add_argument('--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
  args = parser.parse_args()

  if not args.files :
    print "No input files! Nothing to do, so exit."
    quit()

  # Fix negative values of skip and max
  if args.max < 0 :
    args.max = 99999
  if args.skip <= 0 :
    args.skip = -1

  # Read in PDB files
  if len(args.files) == 1 :
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.files[0],skip=args.skip,readmax=args.max)
  else :
    pdbfiles = simulationobjects.PDBSet()
    for filename in args.files[args.skip:args.max+1] :
      pdb = simulationobjects.PDBFile(filename=filename)
      pdbfiles.pdbs.append(pdb)

  print "\nPDB data has been read. Clustering...\n"

  # Find the molecule and cluster.
  confs, numconfs = get_coords(pdbfiles,args.molecule,args.atom)
  molpdb = []
  filenum = 0
  while len(molpdb) == 0:
    molpdb = _GetMolTemplate(args.files[filenum], args.molecule)
    filenum += 1

  if numconfs == 0 and args.atom==None:
    print "\n Molecule '%s' not found in supplied file(s). Please check input. \n" % args.molecule
    quit()
  elif numconfs == 0 and args.atom is not None:
    print "\n Atom '%s' in molecule '%s' not found in supplied file(s). Please check input. \n" % (args.atom,args.molecule)
    quit()

  clusts, occupancies =  cluster_coords(confs,clustmethod=args.type,cutoff=args.cutoff)
  occ_order = np.argsort(occupancies)

  outfile = open(args.out,"w")
  resid = 1
  for i in occ_order[::-1]:
    _printpdb(clusts[i],molpdb,resid,occupancies[i],outfile)
    resid += 1



  


