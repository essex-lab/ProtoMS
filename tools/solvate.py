# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routines to solvate a structure

This module defines a single public function:
solvate
    
Can be executed from the command line as a stand-alone program
"""

import sys
import random

import numpy as np
from scipy.spatial.distance import cdist

import simulationobjects

def solvate(box, ligand=None, protein=None,
            out="solvent_box.pdb", geometry="box",
            padding=10.0, radius=30.0, center="cent",
            namescheme="ProtoMS", log="solvate.log"):

  """Function to solvate ligand and/or protein structures using a
     pre-equilibrated box of waters.

     solvate(box, ligand=None, protein=None,
               out="solvent_box.pdb", geometry="box",
               padding=10.0, radius=30.0, center="cent",
               namescheme="ProtoMS", log="solvate.log")

  Parameters
  ----------
  box : String
      String giving the name of a pdb file containing a box of pre-
      equilibrated water molecules for solvation
  ligand : simulationobjects.PDBFile or string, optional
      Either pdb instance for ligand to be solvated, or string giving
      the name of a pdb file of the ligand. NOTE - either ligand or protein
      or both must be supplied
  protein : simulationobjects.PDBFile or string, optional
      Either pdb instance for protein to be solvated, or string giving
      the name of a pdb file of the protein. NOTE - either ligand or protein
      or both must be supplied
  out : string, optional
      File name to write solvent box. Default of "solvent_box.pdb"
  geometry : string, optional
      Shape of solvent box, options "box", "droplet", "flood". Default box.
      Droplet defines a sphere of radius 'radius' around center 'center'.
      Box defines a box that extends at least distance 'padding' from solute
      Flood defines the equivalent, but waters overlapping with the solute
      are not removed
  padding : float, optional
      Minimum distance between the solute and the box edge. Default 10.0A
  radius : float, optional
      The radius of the added droplet sphere. Default 30.0A
  center : string, optional
      Defines the center of the added droplet sphere. Options "cent",
      one, two or three numbers. Default "cent".
      Cent defines the center based on all solute atom coordinates
      One number defines coordinates for the center as identical x,y,z
      Two numbers defines coordinates for the center over an atom range
      Three numbers difines coordinates for the center as separate x,y,z
  namescheme : string, optional
      Naming scheme for printout, options "ProtoMS", "Amber". Default ProtoMS
  log : string, optional
      Logfile for output, default solvate.log

  Returns
  -------
  None, but writes out a pdb file of the added solvent box,
  including any crystallographic waters in the protein """

  #
  # Lennard-Jones parameters of GAFF atom-types from gaff14.types
  # H and O parameters are averaged over all atom-types
  #
  # This is used to calculate the van der Waals interaction with water oxygen
  #
  params = {}
  params["h"] = (1.71,0.01)
  params["o"] = (3.04,0.19)
  params["c"] = (3.40,0.09)
  params["n"] = (3.25,0.17)
  params["s"] = (3.56,0.25)
  params["p"] = (3.74,0.20)
  params["f"] = (3.12,0.06)
  params["cl"] = (3.47,0.27)
  params["br"] = (3.60,0.31)
  params["i"] = (3.83,0.40)

  # Add buffer to sigmas to avoid over-filling cavities in protein
  # (Equivalent to a 1A increase in vdW radii)
  offset = 0.89

  for element in params:
    params[element] = (params[element][0]+offset,params[element][1])

  #
  # -----------------
  # Utility routines
  # -----------------
  #


  def element(atm) :
    """ Extract the element from a pdb-type atom name
    Parameters
    ----------
       atm - the atom name

    Returns
    -------
       the element
    """
    atm = atm.strip()
    atm = atm.lower()
    if isinstance(atm[0],int) :
      atm = atm[1:]
    if len(atm) > 1 and atm[0:2] == 'cl' : return 'cl'
    if len(atm) > 1 and atm[0:2] == 'om' : return 'om'
    if len(atm) > 1 and atm[0:2] == 'br' : return 'br'
    return atm[0]


  def read_solute(sol_pdb,solute) :
    """ Read a solute pdb-file
    Parameters
    ----------
     sol_pdb - A simulationobjects pdb object of the solute
      solute  - a dictionary of read solute data

    Returns
    -------
     none, but
       solute["atoms"] contains the element of all atoms
       solute["xyz"] contains all the coordinates (Nx3 NumPy array)
   """
    for i in sol_pdb.residues :
      for j in range(len(sol_pdb.residues[i].atoms)) :
        solute["atoms"].append(element(sol_pdb.residues[i].atoms[j].name))
        solute["xyz"].append(sol_pdb.residues[i].atoms[j].coords)

  def read_protein(prot_pdb,solute,solvent) :
    """ Read a protein + xtal waters pdb-file
    Parameters
    ----------
     prot_pdb - a simulationobjects pdb object of the protein + xtal waters
      solute  - a dictionary of read solute data
      solvent  - a dictionary of xtal water data

    Returns
    -------
     none, but
       solute["atoms"] contains the element of all atoms
       solute["xyz"] contains all the coordinates including xtal waters (Nx3 NumPy array)
       solvent["xyz"] contains all the coordinates of just the waters (Nx3 NumPy array)
    """
    solvent["atoms"] = []
    solvent["xyz"] = []
    for i in prot_pdb.residues :
      for j in range(len(prot_pdb.residues[i].atoms)) :
        solute["atoms"].append(element(prot_pdb.residues[i].atoms[j].name))
        solute["xyz"].append(prot_pdb.residues[i].atoms[j].coords)
    for i in prot_pdb.solvents :
      xyz = []
      for j in range(len(prot_pdb.solvents[i].atoms)) :
        solute["atoms"].append(element(prot_pdb.solvents[i].atoms[j].name))
        solute["xyz"].append(prot_pdb.solvents[i].atoms[j].coords)
        xyz.append(prot_pdb.solvents[i].atoms[j].coords)
      logfile.write("Recognised residue %s in protein as a crystallographic water - this is being retained\n" % prot_pdb.solvents[i].index)
      xyz = np.array(xyz)
      solvent["xyz"].append(xyz)

  def read_watbox(filename,watbox) :
    """ Read a pdb-file and extract water coordinates
    Parameters
    ----------
     filename - the name of the file to read
     watbox   - a dictionary of read water data

    Returns
    -------
     none, but
       watbox["xyz"] contains all the coordinates (array of NumPy arrays)
       watbox["min"] is the origin of water box
       watbox["len"] is the length of the water box
    """
    # Take out the water records
    solv_box_pdb = simulationobjects.PDBFile(filename=filename)

    # Then loop over all water molecules
    watbox["xyz"] = []
    minxyz = np.zeros(3)+10000
    maxxyz = np.zeros(3)-10000
    for i in solv_box_pdb.solvents :
      xyz = []
      # And all water atoms
      for j in range(len(solv_box_pdb.solvents[i].atoms)) :
        xyz.append(solv_box_pdb.solvents[i].atoms[j].coords)     
      xyz = np.array(xyz)
      minxyz=np.minimum(minxyz,np.min(xyz,axis=0))
      maxxyz=np.maximum(maxxyz,np.max(xyz,axis=0))
      watbox["xyz"].append(xyz)

    watbox["min"] = minxyz
    watbox["len"] = maxxyz-minxyz   

  #---------------------------------------------------------------- 
  # Routines to remove waters that overlap with solute coordinates
  #----------------------------------------------------------------

  def vdw_energy(element,atomxyz,oxyz) :
    """ Calculate the van der Waals energy between a water oxygen atom and a given atom
    Parameters
    ----------
     element - the element of the atom
     atomxyz - the coordinate of the atom (NumPy array)
     oxyz    - the coordinate of the water oxygen (NumPy array)

    Returns
    -------
     an approximate energy, cut-off at 5 A
    """

    # Return 0 for undefined atoms
    if element not in params : return 0.0

    r2 = np.sum((atomxyz-oxyz)**2)
    # Return 0 beyond a 5 A cut-off
    if r2 > 25.0 : return 0.0

    # Lorent--Berenholz combining rules with a water oxygen
    sigma = 0.5*(params[element][0]+3.15)
    eps = np.sqrt(params[element][1]*0.15)

    sigr2 = sigma*sigma/r2
    return 4.0*eps*(sigr2**6-sigr2**3)


  def remove_overlap(water,solute,cutoff,logfile) :
    """ Remove water molecules that overlap with solute coordinates
    Parameters
    ----------
     water  - the added water molecules (a dictionary)
     solute - the solute molecule (a dictionary)
     cutoff - the energy cutoff to use
     logfile - a file point to write log information

    Returns
    -------
     none, water argument is modified
    """

    logfile.write("\n")
    saved = []
    # Loop over all added water molecules
    for wat in water["xyz"] :
      # Find the smallest distance between the oxygen and any solute atom
      i = np.argmin(cdist(wat,solute["xyz"],'sqeuclidean'),axis=1)[0]
      # Caculate an approximate vdW energy and determine if it should be kept
      ene = vdw_energy(solute["atoms"][i],solute["xyz"][i,:],wat[0,:]) 
      if ene < cutoff :
        saved.append(np.array(wat,copy=True))
      else :
        logfile.write("Overlap: removed water molecule with vdW energy: %.2f\n"%ene)

    logfile.write("Overlap: removed %d water molecules\n"%(len(water["xyz"])-len(saved)))
    water["xyz"] = saved  

  # ---------------------------------------------------------------------
  # Routines to solvate a solute by replicating a pre-equilibrated box
  # ---------------------------------------------------------------------

  def add_template(template,delta,maxbox,container) :
    """ Add a copy of template box to a container of waters
    Parameters
    ----------
  #   template - the template to be copied (a dictionary)
  #   delta    - an displacement to add to the template
  #   maxbox   - the maximum of the big box
  #   container - the container where the template should be added

    Returns
    -------
  #   none, container argument is modified 
    """
    for wat in template :
      wat2 = wat-delta
      if np.all(np.max(wat2,axis=0) <= maxbox) : # Skip water if it is outside box
        container.append(wat2)


  def replicate_box(watbox,solute,padding,flooding,bigbox,logfile) :
    """ Replicate a pre-equilibrated box to cover a solute
    Parameters
    ----------
  #   watbox - the pre-equilibrated box (a dictionary)
  #   solute - the solute (a dictionary)
  #   padding - the minimum distance between solute and box edge
  #   flooding - turn on flooding which defines the bigbox in different way
  #   bigbox  - the added water molecules (a dictionary)
  #   logfile - a file point to write log information

    Returns
    -------
  #  none, bigbox argument is modified
  #
    """
    # Define the extent of the solvation box
    if not flooding :
      bigbox["min"] = solute["cent"] - np.max(solute["len"])/2 - padding
      bigbox["max"] = solute["cent"] + np.max(solute["len"])/2 + padding
    else :
      bigbox["min"] = np.array(solute["min"]) - padding
      bigbox["max"] = np.array(solute["max"]) + padding
    bigbox["xy_plane"] = []
    bigbox["xyz"] = []

    logfile.write("\nReplication: solvation box origin: %.2f %.2f %.2f\n"%(bigbox["min"][0],bigbox["min"][1],bigbox["min"][2]))
    logfile.write("Replication: solvation box length: %.2f %.2f %.2f\n\n"%(bigbox["max"][0]-bigbox["min"][0],bigbox["max"][1]-bigbox["min"][1],bigbox["max"][2]-bigbox["min"][2]))

    # Add boxes in the x-dimension
    addxyz = np.array(bigbox["min"],copy=True)
    while addxyz[0] < bigbox["max"][0] :
      # Add boxes in the y-dimension
      while addxyz[1] < bigbox["max"][1] :
        add_template(watbox["xyz"],watbox["min"]-addxyz,bigbox["max"],bigbox["xy_plane"])    
        logfile.write("Replication: added box at: %7.2f %7.2f %7.2f\n"%(addxyz[0],addxyz[1],addxyz[2]))
        addxyz[1] = addxyz[1] + watbox["len"][1]

      addxyz[1] = bigbox["min"][1]
      addxyz[0] = addxyz[0] + watbox["len"][0]
      logfile.write("\n")

    # Add xy-planes in the z-dimension
    addz = bigbox["min"][2]
    while addz < bigbox["max"][2] :
      delta = np.array([0.0,0.0,bigbox["min"][2]-addz])
      add_template(bigbox["xy_plane"],delta,bigbox["max"],bigbox["xyz"])   
      logfile.write("Replication: added xy-plane at: %7.2f\n"%addz)
      addz = addz + watbox["len"][2]

    logfile.write("\nReplication: added %d water molecules in total\n"%len(bigbox["xyz"]))

  # -------------------------------------------------------------
  # Routines to solvate a solute in a droplet of water molecules 
  # -------------------------------------------------------------
  
  def define_center(cent,solute) :
    """ Parse the center argument given by the user and thereby define the center of the droplet
    Parameters
    ----------
  #   cent - the command-line argument
  #   solute - the solute (NumPy array)

    Returns
    -------
  #   the center of the droplet (NumPy array)
    """
    if cent[0:4] == "cent" :
      return np.array(solute["cent"],copy=True)
    else :
      cols = cent.strip().split()
      if len(cols) == 3 : # Gives the center
        return np.array(cols[:3],dtype=float)
      elif len(cols) == 2 : # Gives an atom range
        # Calculate the center of geometry of atoms between cols[0] and cols[1]
        first = int(cols[0])-1
        last = int(cols[1])
        return np.mean(solute["xyz"][first:last,:],axis=0)
      else : # Gives the center by replication
        return np.array([cols[0],cols[0],cols[0]],dtype=float)

  def rand_sphere(rad) :
    """ Utility to draw a random vector on a sphere
    Parameters
    ----------
  #   rad - the radius of the sphere
  # Output:
  #   the vector (NumPy array)
    """
    x = random.uniform(-rad,rad)
    y = random.uniform(-rad,rad)
    z = random.uniform(-rad,rad)
    while (x**2+y**2+z**2)>rad**2 :
      x = random.uniform(-rad,rad)
      y = random.uniform(-rad,rad)
      z = random.uniform(-rad,rad)  
    return np.array([x,y,z])

  def make_droplet(watbox,solute,radius,droplet,logfile) :
    """ Create a droplet on top of a solute
    Parameters
    ----------
  #   watbox - the pre-equilibrated box (a dictionary)
  #   solute - the solute (a dictionary)
  #   radius - the radius of the droplet
  #   droplet  - the added water molecules (a dictionary)
  #   logfile - a file point to write log information

    Returns
    -------
  #  none, droplet argument is modified
    """
    # Maximum extent if a box covering the entire droplet
    maxxyz = droplet["cent"] + radius
    # Spacing determined from the theoretical density of pure water
    delta  = (1.0/0.0335)**0.3333333

    logfile

    rad2 = radius**2
    droplet["xyz"] = []
    x = droplet["cent"][0] - radius
    while x <= maxxyz[0] :
      y = droplet["cent"][1] - radius
      while y <= maxxyz[1] :
        z = droplet["cent"][2] - radius
        while z <= maxxyz[2] :
          # Check if we are on the sphere
          r2 = (x-droplet["cent"][0])**2+(y-droplet["cent"][1])**2+(z-droplet["cent"][2])**2
          if r2 <= rad2 :
            wi = random.randint(0,len(watbox["xyz"])-1)
            randxyz = rand_sphere(1.0)
            offset = watbox["xyz"][wi][0,:]-np.array([x,y,z])+randxyz
            wat = watbox["xyz"][wi]-offset
            droplet["xyz"].append(wat)
          z = z + delta
        y = y + delta
      x = x + delta

    logfile.write("\nDroplet: radius %.2f\n"%radius)  
    logfile.write("Droplet: center sphere at: %.2f %.2f %.2f\n"%(droplet["cent"][0],droplet["cent"][1],droplet["cent"][2]))  
    logfile.write("Droplet: added %d water molecules in total\n"%len(droplet["xyz"]))

  # -------------------------------------
  # Main routine begins below
  # -------------------------------------

  # Check for presence of required inputs and change pdb files into pdb objects
  if ligand is None and protein is None :
    raise simulationobjects.SetupError("Neither ligand nor protein have been defined to solvate - quitting!")
  if ligand is not None :
    if isinstance(ligand,simulationobjects.PDBFile) :
      pass
    else :
      ligand = simulationobjects.PDBFile(filename=ligand)
  if protein is not None :
    if isinstance(protein,simulationobjects.PDBFile) :
      pass
    else :
      protein = simulationobjects.PDBFile(filename=protein)


  # Define water naming scheme
  names = []
  resname = []
  if namescheme == "ProtoMS" :
    names = ["O00 ","H01 ","H02 ","M03 "]
    resname = ["","","T3P","T4P"]
  elif namescheme == "Amber" :
    names = [" O  "," H1 "," H2 ","EPW "]
    resname = ["","","WAT","WAT"]

  # Open logfile
  logfile = open(log,"w")

  # Read the pre-equilibrated water box
  watbox = {}
  read_watbox(box,watbox)

  # Read the solute coordinates
  solute = {}
  solute["atoms"] = []
  solute["xyz"] = []
  solvent = {}
  if ligand is not None :
    read_solute(ligand,solute)
  if protein is not None :
    read_protein(protein,solute,solvent)
    solvent["xyz"] = np.array(solvent["xyz"])
  solute["xyz"] = np.array(solute["xyz"])
  solute["min"] = np.min(solute["xyz"],axis=0)
  solute["max"] = np.max(solute["xyz"],axis=0)
  solute["len"] = solute["max"] - solute["min"]
  solute["cent"] = np.average(solute["xyz"],axis=0)

  # Write some log information
  logfile.write("Water box min %.2f %.2f %.2f\n"%(watbox["min"][0],watbox["min"][1],watbox["min"][2]))
  logfile.write("Water box len %.2f %.2f %.2f\n"%(watbox["len"][0],watbox["len"][1],watbox["len"][2]))
  logfile.write("Solute min %.2f %.2f %.2f\n"%(solute["min"][0],solute["min"][1],solute["min"][2]))
  logfile.write("Solute max %.2f %.2f %.2f\n"%(solute["max"][0],solute["max"][1],solute["max"][2]))
  logfile.write("Solute cent %.2f %.2f %.2f\n"%(solute["cent"][0],solute["cent"][1],solute["cent"][2]))

  # Now add waters
  added_water = {}
  if geometry == "box" :
    replicate_box(watbox,solute,padding,False,added_water,logfile)
  elif geometry == "droplet" :  
    added_water["cent"] = define_center(center,solute)
    make_droplet(watbox,solute,radius,added_water,logfile)
  elif geometry == "flood" :
    replicate_box(watbox,solute,0.1,True,added_water,logfile)

  # Remove water that overlap with the solute
  if not geometry == "flood" :
    remove_overlap(added_water,solute,20.0,logfile)
  new_watbox = simulationobjects.PDBFile()
  
  # Write box or cap information as header
  if geometry == "box" :
    header = "HEADER box %.4f %.4f %.4f %.4f %.4f %.4f\n"%(added_water["min"][0]-0.5,added_water["min"][1]-0.5,added_water["min"][2]-0.5,added_water["max"][0]+0.5,added_water["max"][1]+0.5,added_water["max"][2]+0.5)
  elif geometry == "droplet" :
    header = "HEADER cap %.4f %.4f %.4f %.4f 1.5\n"%(added_water["cent"][0],added_water["cent"][1],added_water["cent"][2],radius)
  try:
    for i in range(len(solvent["xyz"])) :
      if len(solvent["xyz"][i]) == len(added_water["xyz"][0]) :
        added_water["xyz"].append(solvent["xyz"][i])
      else :
        raise simulationobjects.SetupError("The length of a crystal water molecule does not match the length of molecules in the water box - please check they are the same type (TIP3P/TIP4P)!")
        quit()
  except KeyError:
    logfile.write("No crystallographic waters added to solvent box!\n")
    pass

  # Write water coordinates
  atmidx = 1
  for i,w in enumerate(added_water["xyz"]) :
    resid = i+1
    if resid >= 10000 : resid = resid - 9999
    if atmidx >= 100000 : atmidx = atmidx - 99999
    new_watbox.solvents[resid] = simulationobjects.Residue(name=resname[len(w)-1],index=resid)
    newatom1 = simulationobjects.Atom(index=atmidx,name=names[0],resindex=resid,resname=resname[len(w)-1],coords=w[0])
    atmidx += 1
    new_watbox.solvents[resid].addAtom(atom=newatom1)
    newatom2 = simulationobjects.Atom(index=atmidx,name=names[1],resindex=resid,resname=resname[len(w)-1],coords=w[1])
    atmidx += 1
    new_watbox.solvents[resid].addAtom(atom=newatom2)
    newatom3 = simulationobjects.Atom(index=atmidx,name=names[2],resindex=resid,resname=resname[len(w)-1],coords=w[2])
    atmidx += 1
    new_watbox.solvents[resid].addAtom(atom=newatom3)
    if len(w) > 3 :
      newatom4 = simulationobjects.Atom(index=atmidx,name=names[3],resindex=resid,resname=resname[len(w)-1],coords=w[3])
      atmidx += 1
      new_watbox.solvents[resid].addAtom(atom=newatom4)

  new_watbox.write(filename=out, header=header)
  logfile.close()



#
# -------------------------------------
# If this is run from the command-line
# -------------------------------------
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  disclaimer = """
    if -b or -s are not supplied on the command-line, the program will ask for them.    \n
    -c can be either 'cent' or a string containing 1, 2 or 3 numbers. If 1 number is
    given it will be used as center of the droplet in x, y, and z. If 2 numbers are 
    given this is interpreted as an atom range, such that the droplet will be centered 
    on the indicated atoms, and if 3 numbers are given this is directly taken as the
    center of droplet  \n     
    Example usages:
      solvate.py -b ${PROTOMSHOME}/tools/sbox1.pdb -s solute.pdb
         (will solvate 'solute.pdb' in a box that extends at least 10 A from the solute)
      solvate.py -b ${PROTOMSHOME}/tools/sbox1.pdb -s protein.pdb -g droplet -r 25.0
         (will solvate 'protein.pdb' in a 25 A droplet centered on all coordinates)
  """
  parser = argparse.ArgumentParser(description="Program to solvate a solute molecule in either a box or a droplet",epilog=disclaimer,formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-b','--box',help="a PDB-file containing a pre-equilibrated box of water molcules",default="")
  parser.add_argument('-s','--solute',help="a PDB-file containing the solute molecule",default=None)
  parser.add_argument('-pr','--protein',help="a PDB-file containing the protein molecule",default=None)
  parser.add_argument('-o','--out',help="the name of the output PDB-file containing the added water, default solvent_box.pdb",default="solvent_box.pdb")
  parser.add_argument('-g','--geometry',choices=["box","droplet","flood"],help="the geometry of the added water, should be either 'box', 'droplet' or 'flood'",default="box")
  parser.add_argument('-p','--padding',type=float,help="the minimum distance between the solute and the box edge, default=10 A",default=10.0)
  parser.add_argument('-r','--radius',type=float,help="the radius of the droplet, default=30A",default=30.0)
  parser.add_argument('-c','--center',help="definition of center, default='cent'",default="cent")
  parser.add_argument('-n','--names',choices=["Amber","ProtoMS"],help="the naming convention, should be either Amber or ProtoMS",default="ProtoMS")
  parser.add_argument('-l','--log',help="The logfile for output, default is solvate.log",default="solvate.log")
  args = parser.parse_args()

  # Ask for input that is absolutely necessary and that do not have any defaults
  if args.solute is None and args.protein is None:
    print "You haven't entered a protein or a solute to solvate!\n"
    print "Please enter a solute now (Return to skip):"
    args.solute = raw_input()
    print "Please enter a protein now (Return to skip):"
    args.protein = raw_input()
    if args.solute == "" and args.protein == "":
      print "You still haven't entered a protein or a solute to solvate!\n"
      quit()
  if args.box == "" :
    print "Enter the filename of a pre-equilibrated water box: ",
    args.box = raw_input()
  if args.solute == "" : args.solute = None
  if args.protein == "" : args.protein = None

  solvate(args.box,args.solute,args.protein,args.out,args.geometry,args.padding,args.radius,args.center,args.names,args.log)

