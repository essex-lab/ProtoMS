# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Julien Michel
#          Gregory Ross

""" 
  Useful classes for setup and analysis of ProtoMS simulations.
  Contain classes to 
    1) Read, modify and write structures in PDB format
    2) Store parameter collections
    3) Read and store ProtoMS results files
    4) Read, modify and write ProtoMS template files
"""

import math
import sys
import os
import glob
import logging
import copy
import gettext
import argparse

import numpy as np

boltz = 0.00198717076 # kcal.mol-1K-1

def is_solvent(name) :
  """
  Check whether a given residue name
  is a solvent residue name
  
  Parameters
  ----------
  name : string
    the residue name
 
  Return
  ------
  Boolean
    whether the residue name
    is a solvent residue name
  """
  return name in ['WAT','wat','HOH','hoh','DOD','dod','T3P','t3p','T4P','t4p','SOL','sol','see','SEE']
    

class SetupError(Exception) :
    """ A general exception class to be raised by setup code
    """
    pass

#------------------------------------------
# Classes and routines to handle PDB-files
#------------------------------------------

class Atom(object):
    """ 
    Class for holding a PDB atom record

    Attributes
    ----------
    index : integer
      the number of the atom in the pdb file
    resindex : integer
      the number of the residue in the pdb file
    name : string
      the name of the atom
    resname : string
      the name of the residue
    coords : NumpyArray
      Cartesian coordinates
    type : string
      either atom or hetatm
    element : string
      the element of the atom
    """
    def __init__(self,index=0,name="?",resname="???",resindex=0,coords=[]):
        self.index = index
        self.resindex = resindex
        self.name = name
        self.resname =resname
        self.coords = np.array ( coords )
        self.type = "atom"
        self.element = "??"      
    def getElement(self):
        """ 
        Set the element of this atom        
        """
        k = 0
        name = self.name.strip()
        self.element = name[k]
        while self.element.isdigit():
            k = k + 1
            self.element = name[k]
        if len(name) > k+1 :
            if name[k].lower() + name[k+1].lower() in ["cl","br","mg","mn"]:
                self.element = name[k].lower() + name[k+1].lower()
        else:
            self.element = name[k].lower()
    def __str__(self):
        """
        Produces a string representation, viz. a standard ATOM record
        """
        return "ATOM  %-5d %3s %3s    %4d     %8.3f%8.3f%8.3f  1.00 0.00" % (self.index,self.name,self.resname,self.resindex,self.coords[0],self.coords[1],self.coords[2])
              
class Residue(object):
    """ 
    Class for holding a set of PDB atoms

    Attributes
    ----------
    atoms : list of Atom objects 
      the atom of this residue
    index : integer
      the number of the residue in the pdb file
    name : string
      the name of the residue
    center : NumPy array
      the center of coordinates of all atoms
    """
    def __init__(self,name="???",index=0):
        self.atoms = []
        self.index = index
        self.name = name
        self.center = np.array([0.0,0.0,0.0])          
    def addAtom(self,atom):
        """ 
        Adds an atom to this residue

        Parameters
        ----------
        atom : Atom object
          the atom to add
        """        
        assert atom.type ==  "atom"
        self.atoms.append(atom)
    def getCenter ( self ):
        """ 
        Sets the center of coordinates for this residue 
        """
        coords = np.array ( [ atom.coords for atom in self.atoms ] )
        self.center = coords[:,0].mean(), coords[:,1].mean(), coords[:,2].mean()
        #return self.center            
    def isAminoacid (self):
        """
        Checks is the residue name is an aminoacid name
        """
        return self.name.upper() in  ['GLH', 'ILE', 'GLN', 'GLY', 'GLU', 'HIP', 'HIS', 'SER', 'LYS', 'PRO', 'CYX', 'HIE', 'LYN', 'ASH', 'ASN', 'CYS', 'VAL', 'THR', 'ASP', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
    def __str__(self):
        """
        Produces a string representation, viz. lines of ATOM records
        """
        return '\n'.join(atom.__str__() for atom in self.atoms)

class PDBFile:
    """ 
    Class for holding the atoms and residues of a PDB file
     
    Attributes
    ----------
    center : NumPy array
      the center of coordinates of all atoms
    header : string
      additional lines to print before ATOM records
    name : string
      the name of the file
    residues : dictionary of Residue objects
      all non-solvent residues
    solvents : dictionary of Residue objects
      all solvent residues
    """
    def __init__ ( self, filename = None ):
        self.residues = {}
        self.solvents = {}
        self.name = ""
        self.center = np.array([0.0,0.0,0.0])
        self.header = ""

        if filename is not None : self.read(filename)
    def __str__ ( self ):
        """
        Returns the filename
        """
        return self.name
    def copy ( self ) :
        """ 
        Make a copy of the residues and solvents dictionaries, 
        not the Residue and Atom objects themselves

        Returns
        -------
        PDBFile
          the newly created instance
        """
        new = PDBFile()
        new.name = self.name
        new.header = self.header
        new.center = np.array(self.center,copy=True)
        for res in self.residues : new.residues[res] = self.residues[res]
        for sol in self.solvents : new.solvents[sol] = self.solvents[sol]
        return new
    def read(self, filename) :
        """
        Read a pdb-file from disc
       
        It will overwrite existing residues and renumber all residues from 1
       
        Parameters
        ----------
        filename : string
          the name of the pdb file
        """
        self.name = filename
        with open ( filename ) as f:
          self.read_from(f)
    def read_from(self,f,resname=None) :
        """
        Read a pdb structure from a fileobject
        
        Terminates reading when found a record starting with END
        
        Parameters
        ----------
        f : file object
          the object to read from
        resname : string, optional
          the name of the residue to read, only read this        
        """
        residues = {}
        solvents = {}
        line = f.readline()
        nres = 0
        prevres = -1
        while line :
          if line[:6] in ["ATOM  ","HETATM"] :           
              restype = line[17:20]
              # Check if should skip this line
              if resname is not None and restype.lower() != resname.lower() : 
                line = f.readline()
                continue

              index = int(line[6:11].strip())
              atname = line[12:16].strip()
              resnum = int(line[22:26].strip())
              if resnum != prevres :
                nres = nres + 1
              prevres = resnum
              x,y,z = float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())
              coords = [x,y,z]
              newatom = Atom(index=index,name=atname,resindex=nres,
                             resname=restype,coords=coords)
              # If solvent
              if not is_solvent(restype) :
                  try:
                      residues[nres]
                  except KeyError:
                      residues[nres] = Residue(name=restype,index=nres)
                  #print "resnum is %d adding atom %d" % (resnum,index)
                  residues[nres].addAtom(atom=newatom)
              else:
                  #print restype
                  try:
                      solvents[nres]
                  except KeyError:
                      solvents[nres] = Residue(name=restype,index=nres)
                  if len (solvents[nres].atoms) >= 4:
                      raise BaseException("More than one residue with number %d " % nres)
                  solvents[nres].addAtom(atom=newatom)
                  #print solvents
          elif line[:3] == "TER" :
            nres = nres + 1
          elif line[:6] in ["HEADER","REMARK"] :
            self.header = self.header + line
          elif line[:3] == "END" :
            self.residues,self.solvents = residues,solvents
            return
          # Read next line
          line = f.readline()
        self.residues,self.solvents = residues,solvents
    def write ( self, filename, renumber = False, header = None, solvents = True ):
        """
        Write the PDB file to disc

        Parameters
        ----------
        filename : string
          the name of the file
        renumber : boolean, optional
          indicate if residues should be renumbered
        header : string, optional
          additional lines to print before the atom records
        solvents : boolean, optional
          if to write the solvents
        """
        with open ( filename, 'w' ) as f:
            if header is None :
                f.write ( self.header)
            else :
                f.write ( header )
            for i, res in enumerate ( sorted ( self.residues.keys() ), 1 ):
                for atom in self.residues[res].atoms:
                    if renumber:
                        atom.resindex = i
                    s = "ATOM  %5d %-4s %3s  %4d    %8.3f%8.3f%8.3f        \n" % (atom.index,atom.name,atom.resname,atom.resindex,atom.coords[0],atom.coords[1],atom.coords[2])
                    f.write ( s )
                if i < len(self.residues.keys()) and not self.residues[self.residues.keys()[i]].isAminoacid() : f.write ( "TER \n" )     
            if len(self.residues.keys()) > 0 and len(self.solvents.keys()) > 0 and solvents : f.write ( "TER \n" )
            if solvents :
                l = sorted ( self.solvents.keys() )
                for i, sol in enumerate (l , 1 ):
                    for atom in self.solvents[sol].atoms:
                        if renumber:
                            atom.resindex = i
                        s = "ATOM  %5d %-4s %3s  %4d    %8.3f%8.3f%8.3f        \n" % (atom.index+1,atom.name,atom.resname,atom.resindex,atom.coords[0],atom.coords[1],atom.coords[2])
                        f.write ( s )
                    if i < len(l) : f.write ( "TER \n" )     
    def getCenter(self):
        """
        Calculate the center of geometry  
        
        Returns
        -------
        NumPy array
          the center
        """
        center = np.array([0.0,0.0,0.0])
        for res in self.residues:
            self.residues[res].getCenter()
            rescent = self.residues[res].center
            center = center + rescent
        self.center = center/float(len(self.residues))
        return self.center
    def getBox(self,atomlist = 'all',reslist = 'all') :
        """
        Calculate the smallest box to encompass the atoms

        Parameters
        ----------
        atomlist = list, optional
          the atoms to be taken into acount to get the box
        reslist = list, optional
          the residues to be taken into acount to get the box

        Returns
        -------
        dictionary of Numpy arrays
          the center, length and origin of the box
        """
        if reslist is 'all' :
          reslist = list(set([self.residues[i].name for i in self.residues] + [self.solvents[i].name for i in self.solvents]))	 
        if atomlist is 'all' :
          atomlist = list(set([atom.name for i in self.residues for atom in self.residues[i].atoms] + [atom.name for i in self.solvents for atom in self.solvents[i].atoms]))
        minxyz = np.zeros(3)+1E6
        maxxyz = np.zeros(3)-1E6
        for res in self.residues :
          for atom in self.residues[res].atoms :
            if atom.name in atomlist and self.residues[res].name in reslist :
              minxyz = np.minimum(minxyz,atom.coords)
              maxxyz = np.maximum(maxxyz,atom.coords)
        for res in self.solvents :
          for atom in self.solvents[res].atoms :
            if atom.name in atomlist and self.solvents[res].name in reslist :
              minxyz = np.minimum(minxyz,atom.coords)
              maxxyz = np.maximum(maxxyz,atom.coords)
        return {"center":minxyz+(maxxyz-minxyz)/2.0,"len":maxxyz-minxyz,"origin":minxyz}
       
class  PDBSet :
    """
    Hold a collection of PDBFile objects
    
    Attributes
    ----------
    pdbs : list of PDBFile objects
      the pdb files
    """
    def __init__(self) :
        self.pdbs = []
    def from_residues(self,pdbfile) :
        """
        Split a PDB file into a set by residue
        
        Parameters
        ----------
        pdbfile : PDBFile object
          the pdb-file to split          
        """
        self.pdbs = []
        for resi,res in pdbfile.residues.iteritems() :
          self.pdbs.append(PDBFile()) 
          self.pdbs[-1].residues[resi] = res
          self.pdbs[-1].header = pdbfile.header
        for soli,sol in pdbfile.solvents.iteritems() :
          self.pdbs.append(PDBFile()) 
          self.pdbs[-1].solvents[soli] = sol
          self.pdbs[-1].header = pdbfile.header
    def read(self,filename,resname=None) :
        """
        Read a set of pdb structures from a file
        
        Parameters
        ----------
        filename : string
          the name of the file to read
        resname : string, optional
          the name of the residue to read, only read this
        """
        self.pdbs = []
        with open(filename,"r") as f :
            while True :
                pdb = PDBFile()
                pdb.read_from(f,resname=resname)
                if not (pdb.residues or pdb.solvents) :
                    break
                else :
                    self.pdbs.append(pdb)
    def write(self,filenames,solvents=True) :
        """
        Write the set to disc

        If one filename is given all PDB-files are written to one file
        
        Parameters
        ----------
        filenames : list of strings
          the name of the file(s) to write to
        solvents : boolean
          if to write the solvents
        """        

        if len(filenames) == len(self.pdbs) :
          for pdb,filename in zip(self.pdbs,filenames) :
            pdb.write(filename,solvents=solvents)
        elif isinstance(filenames,basestring) or len(filenames) == 1 :
          raise SetupError("This feature is not implemented yet")
        else :
          raise SetupError("Invalid number of filenames given")        
  
def merge_pdbs(pdbobjs) :
    """
    Merge pdb files

    Create a new PDBFile instance and add all residues and solvents from
    all pdb files given, renumbering residues

    Parameters
    ----------
    pdbobjs : list of PDBFile objects
 
    Returns
    -------
    PDBFile object
      the newly created and merged pdb structure
    """
    pdbout = PDBFile()
    nres = 0
    names = []
    for pdbobj in pdbobjs :
        names.append(pdbobj.name)
        for res in sorted ( pdbobj.residues.keys() ) :
          nres = nres + 1
          pdbobj.residues[res].index = nres
          for atom in pdbobj.residues[res].atoms : atom.resindex = nres
          pdbout.residues[nres] = pdbobj.residues[res]
        for sol in sorted ( pdbobj.solvents.keys() ) :
          nres = nres + 1
          pdbobj.solvents[sol].index = nres
          for atom in pdbobj.solvents[sol].atoms : atom.resindex = nres
          pdbout.solvents[nres] = pdbobj.solvents[sol]
    pdbout.name = "_".join(names)
    
    return pdbout      

def find_box(pdbobj) :
  """
  Parse box information from a pdb file

  Looks for CENTER, DIMENSIONS and box keywords
  in the header of the pdb file

  Parameters
  ----------
  pdbobj : PDBFile object
    the structure to parse

  Returns
  -------
  dictionary of Numpy arrays
    the center, length and origin of the box  
  """
  if pdbobj.header == "" : return None
  center = dim = box = None
  headerlines = pdbobj.header.split("\n")
  for line in headerlines :
    if line.find("CENTER") > -1 :
      center = np.array(line.strip().split()[5:],float)            
    elif line.find("DIMENSIONS") > -1 :
      dim = np.array(line.strip().split()[5:],float)
    elif line.find("HEADER box") > -1  or line.find("REMARK box") > -1 :
      box = np.array(line.strip().split()[2:],float)
  if center is not None and dim is not None :
    return {"center":center,"len":dim}
  elif box is not None :
    return {"origin":box[:3],"len":box[3:]-box[:3]}
  else :
    return None
    
def write_box(filename,box) :
  """
  Write a box in PDB file format
  
  Parameters
  ----------
  filename : string
    the name of the file to write
  box : dictionary of Numpy array
    the box specification
  """
  def makeBoundary(box_min,box_max) :
    c1 = (box_min[0],box_min[1],box_min[2]) # Origin
    c2 = (box_max[0],box_min[1],box_min[2])
    c3 = (box_max[0],box_min[1],box_max[2]) 
    c4 = (box_min[0],box_min[1],box_max[2])
    c5 = (box_min[0],box_max[1],box_min[2])
    c6 = (box_max[0],box_max[1],box_min[2])
    c7 = (box_max[0],box_max[1],box_max[2]) 
    c8 = (box_min[0],box_max[1],box_max[2])

    return (c1,c2,c3,c4,c5,c6,c7,c8)

  if "center" not in box :
    box["center"] = box["origin"] + box["len"]/2.0
  elif "origin" not in box :
    box["origin"] = box["center"] - box["len"]/2.0
  boundary = makeBoundary(box["origin"],box["origin"]+box["len"])

  with open(filename,'w') as f :
    f.write("HEADER    CORNERS OF BOX\n")
    f.write("REMARK    CENTER (X Y Z)   %.3f  %.3f  %.3f\n"%(box["center"][0],box["center"][1],box["center"][2]))
    f.write("REMARK    DIMENSIONS (X Y Z)   %.3f  %.3f  %.3f\n"%(box["len"][0],box["len"][1],box["len"][2]))
    f.write("ATOM      1  DUA BOX     1    %8.3f%8.3f%8.3f\n"%boundary[0])
    f.write("ATOM      2  DUB BOX     1    %8.3f%8.3f%8.3f\n"%boundary[1])
    f.write("ATOM      3  DUC BOX     1    %8.3f%8.3f%8.3f\n"%boundary[2])
    f.write("ATOM      4  DUD BOX     1    %8.3f%8.3f%8.3f\n"%boundary[3])
    f.write("ATOM      5  DUE BOX     1    %8.3f%8.3f%8.3f\n"%boundary[4])
    f.write("ATOM      6  DUF BOX     1    %8.3f%8.3f%8.3f\n"%boundary[5])
    f.write("ATOM      7  DUG BOX     1    %8.3f%8.3f%8.3f\n"%boundary[6])
    f.write("ATOM      8  DUH BOX     1    %8.3f%8.3f%8.3f\n"%boundary[7])
    f.write("CONECT    1    2    4    5\n")
    f.write("CONECT    2    1    3    6\n")
    f.write("CONECT    3    2    4    7\n")
    f.write("CONECT    4    1    3    8\n")
    f.write("CONECT    5    1    6    8\n")
    f.write("CONECT    6    2    5    7\n")
    f.write("CONECT    7    3    6    8\n")
    f.write("CONECT    8    4    5    7\n")
    
#--------------------------------
# Classes to hold parameter sets
#--------------------------------

class Parameter:
    """ 
    Class to hold a parameter from a ProtoMS template file
    
    Not valid for dihedral parameters but contains sufficient
    information to check that a parameter exists

    Attributes
    ----------
    index : integer
      the serial number of the parameter
    ats : list of string
      the name of the atoms associated with the parameter
    k : float
      the force constant
    b0 : float
      the equilibrium value
    """
    def __init__ ( self, index, ats, k, b0 ):
        self.index = int ( index )
        self.ats = ats
        self.k = float ( k )
        self.b0 = float ( b0 )

class ParameterSet:
    """
    Class to hold a collection of parameters
    
    Attributes
    ----------
    file_name : string
      the name of the file where the parameters originated
    params : list of Parameter objects
      the parameters in this collection
    """
    def __init__ ( self, param_file, ptype ):
        self.file_name = param_file
        
        assert ptype in ['bond','angle','dihedral']
        pars = []
        atms = []
        with open ( param_file ) as f:
            while not f.next().startswith ( 'mode %s' % ptype ):
                pass
            
            for line in f:
                if line.startswith ( "par" ):
                    pars += [ line ]
                elif line.startswith ( "atm" ):
                    atms += [ line ]
                elif line.startswith ( "mode" ):
                    break
        
        self.params = []
        for par, at in zip ( pars, atms ):
            par_cols, at_cols = par.split(), at.split()
            if par_cols[1] != at_cols[-1]:
                raise ValueError ( "No matching values in gaff parameter file" )
            self.params += [ Parameter ( par_cols[1], at_cols[1:-1], par_cols[2], par_cols[3] ) ]

    def get_params ( self, ats ):
        """
        Find parameters for specific atoms
        
        Parameters
        ----------
        ats : list of string
          the query atoms

        Returns
        -------
        Parameter object
          the found parameter
        """
        try:
            return [ i for i in self.params 
                     if i.ats == ats or i.ats == ats[::-1] ][0]
        except IndexError:
            return [ i for i in self.params 
                     if i.ats[1:3] in [ ats[1:3], ats[1:3:-1] ] ][0]
#--------------------------------------------------
# Classes to read and store ProtoMS restart files
#--------------------------------------------------

class RestartFile :
  """
  Class to hold the information on a
  restart file
  
  Attributes
  ----------
  filename : string
    the name of the restart file
  nsolutes : int
    the number of solute molecules
  ngcsolutes : int
    the number of gcsolute molecules
  solutesdic : dict
    the correlation of solute ids
    with residue names
  gcsolutesdic : dict
    the correlation of gcsolutes ids
    with residue names
  """
  def __init__(self,filename=None) : 
    self.filename = ""
    self.nsolutes = None
    self.ngcsolutes = None
    self.solutesdic = {}
    self.gcsolutesdic = {}
    if filename is not None :
      self.read(filename=filename)

  def read(self,filename) :
    """
    Read the restart file from disc
    
    Parameters
    ----------
    filename : string
      the name of the file to read
    """

    with open(filename,"r") as f :
      for line in f.readlines() :
        cols = line.strip().split()
        if 'nsolutes' in cols[0].lower() :
          self.nsolutes = int(cols[1])
        elif 'gcsolutes' in cols[0].lower() :
          self.ngcsolutes = int(cols[1])
        elif 'gc-solute' in cols[0].lower() :
          self.gcsolutesdic[int(cols[1])] = cols[3]
        elif 'solute' in cols[0].lower() :
          self.solutesdic[int(cols[1])] = cols[3]


#--------------------------------------------------
# Classes to read and store ProtoMS results files
#--------------------------------------------------

class EnergyResults :
  """
  Class to hold energy results
  
  Attributes
  ----------
  curr : float
    the energy at the current lambda
  back : float
    the energy at the previous lambda
  forward : float
    the energy at the next lambda
  type : string
    a label
  """
  def __init__(self,line=None) : 
    self.curr = 0.0
    self.back = 0.0
    self.forw = 0.0
    self.type = ""
    if line is not None :
      self.parse_line(line)
  def parse_line(self,line) :
    """
    Parse energies from the result file
    
    Parameters
    ----------
    line : string
      the line to parse
    """
    cols = line.strip().split()
    self.type = cols[0] 
    self.curr = float(cols[1])
    self.forw = float(cols[6])
    self.back = float(cols[10])
  def __str__(self) :
    if isinstance(self.curr,float) :
      return "%s %20.10F %20.10F %20.10F"%(self.type,self.curr,self.back,self.forw)
    else :
      return "%s Numpy array with %d elements"%(self.type,self.curr.shape[0])
  def __add__(self,b) :
    r = EnergyResults()
    r.curr = self.curr + b.curr
    r.forw = self.forw + b.forw
    r.back = self.back + b.back
    r.type = self.type + b.type
    return r
  def __radd__(self,b) :
    r = EnergyResults()
    r.curr = self.curr + b.curr
    r.forw = self.forw + b.forw
    r.back = self.back + b.back
    r.type = b.type + self.type
    return r

class SnapshotResults :
  """
  Class to store a single snapshot of results
  
  Not all attributes might not be set as they might not
  exists in the result file

  internal_energies and interaction_energies are dictionary
  where the key is is the molecules involved in the energy, e.g.
  protein1, or protein1-solute1. The value is a list of EnergyResults
  objects, for the various energy components.

  Attributes
  ----------
  lam : float
    the current lambda value
  lamb : float
    the previous lambda value
  lamf : float
    the next lambda value
  datastep : int
    the number of steps the averages are calculate over
  lambdareplica : int
    the index of the lambda replica
  temperature : float
    the temperature of the simulation
  ngcsolutes : int
    the number of GC solutes
  nthetasolutes : int
    the number of theta solutes
  bvalue : float
    the Adams value
  solventson : int
    the number of GC solutes that are turned on
  pressure : float
    the pressure
  volume : float
    the volume
  seed : int
    the random seed
  backfe : float
    the free energy to the previous lambda
  forwfe : float
    the free energy to the next lambda
  total : EnergyResults object
    the total energy
  internal_energies : dictionary of lists of EnergyResults object
    the internal energies
  interaction_energies : dictionary of lists of EnergyResults objects
    the interaction energies
  capenergy : EnergyResults object
    the energy of the cap
  extraenergy : EnergyResults object
    all extra energies
  feenergies : dictionary of float
    the total energy at various lambda values
  gradient : float
    the numerical gradient of lambda
  agradient : float
    the analytical gradient of lambda
  thetavals : list of float
    the theta values of all GC solutes
  thetasolvals : list of float
    the theta values of all theta solutes
  """
  def __init__(self) :
    self.thetavals = []
    self.thetasolvals = []
  def parse(self,fileobj) : 
    """
    Parse a file for results
    
    Parameters
    ----------
    fileobj : file object
      the file to read from

    Returns
    -------
    string :
      the line last read
    """
    # First look for the lambda values
    line = fileobj.readline()
    while line.find("RESULTS for lambda") == -1 : line = fileobj.readline()
    cols = line.strip().split()
    if len(cols) == 5 :
      self.lam = float(cols[4])
    else :
      self.lamb = float(cols[3])
      self.lam = float(cols[5])
      self.lamf = float(cols[7])
    
    # Extract various energies and properties
    line = fileobj.readline()
    line = fileobj.readline()
    while line[0] != "#" :
      if line.startswith(" Number of data steps") : self.datastep = int(line.split("=")[1].strip()) 
      if line.startswith(" Lambda replica") : self.lambdareplica = int(line.split("=")[1].strip())
      if line.startswith(" Temperature") : self.temperature = float(line.split("=")[1].strip().split()[0])
      if line.startswith(" Solvents,Proteins,GC-solutes") : self.ngcsolutes =  int(line.split("=")[1].strip().split()[2])
      if line.startswith(" Simulation B factor") : self.bvalue = float(line.split("=")[1].strip())
      if line.startswith(" Simulation B value") : self.bvalue = float(line.split("=")[1].strip())
      if line.startswith(" Molecules in grid") : self.solventson = int(line.split("=")[1].strip())
      if line.startswith(" Pressure") : self.pressure = float(line.split("=")[1].strip().split()[0])
      if line.startswith(" Volume") : self.volume = float(line.split("=")[1].strip().split()[0])
      if line.startswith(" Random number seed") : self.seed = int(line.split("=")[1].strip())
      if line.startswith("Backwards Free Energy") : self.backfe = float(line.split("=")[1].strip().split()[0])
      if line.startswith("Forwards  Free Energy") : self.forwfe = float(line.split("=")[1].strip().split()[0])
      line = fileobj.readline()

    # Look for the total energy
    while not line.startswith("Total Energy ") : line = fileobj.readline()
    self.total = EnergyResults(line=line.replace("Total Energy","Total"))
    
    # Look for internal and average energies
    
    self.internal_energies = {}
    self.interaction_energies = {}
    while line :
      if line.startswith("Internal ") or  line.startswith(" Internal ") :
        cols = line.strip().split()  
        if cols[4] == "solute" :
          key = cols[6]
        else :
          key = cols[4]+cols[5]
        self.internal_energies[key] = []
        line = fileobj.readline() # Dummy line
        line=fileobj.readline()
        while line[0] != "#" :
          self.internal_energies[key].append(EnergyResults(line=line)) # Coul
          line = fileobj.readline() # Dummy line    
          line=fileobj.readline()
      elif line.startswith("Average extra energies") :
        pass
      elif line.startswith(" Average inter-molecular interaction energies") :
        pass
      elif line.startswith("Average total extra energy") :
        line = fileobj.readline() # Dummy line
        self.extraenergy = EnergyResults(line=fileobj.readline())
      elif line.startswith("Average solvent cap energy") :
        self.capenergy = float(line.strip().split()[4])
      elif line.startswith("Average ") :
        cols = line.strip().split()
        if cols[1] in ["solvent-solvent","GCS-GCS"] :
          key = cols[1] 
        elif cols[1] == "protein-protein" :
          key = "protein"+cols[4]+"-protein"+cols[7]
        elif cols[1] == "solute-protein" :
          key = "protein"+cols[4]+"-"+cols[8]
        elif cols[1] == "protein-solvent" :
          key = "protein"+cols[4]+"-solvent"
        elif cols[1] == "solute-solute" :
          key = cols[5]+"-"+cols[8]
        elif cols[1] == "solute-solvent" :
          key = cols[5]+"-solvent"
        elif cols[1] == "solute-GCS" :
          key = cols[5]+"-GCS"   
        self.interaction_energies[key] = []   
        line = fileobj.readline() # Dummy line
        self.interaction_energies[key].append(EnergyResults(line=fileobj.readline())) # Coul
        line = fileobj.readline() # Dummy line    
        self.interaction_energies[key].append(EnergyResults(line=fileobj.readline())) # LJ
      elif line.startswith("FREE ENERGY DATA") :
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # Lambda header
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # First lambda
        self.feenergies = {}
        while line[0] != "#" :
          cols = line.strip().split()
          if cols[2] == "************" :
            cols[2] = "inf"
          self.feenergies[float(cols[0])] = float(cols[2])
          line = fileobj.readline()
        cols = fileobj.readline().strip().split()  
        self.gradient = float(cols[1])
        if len(cols) > 2 :
          self.agradient = float(cols[2])
      elif line.startswith(" Individual theta values"):
        self.thetavals = [None]*self.ngcsolutes
        line = fileobj.readline() # Dummy line
        for i in range(self.ngcsolutes) :
          self.thetavals[i] = (fileobj.readline().strip().split()[2])
      elif line.startswith(" Individual theta solute values"):
        self.nthetasolutes = int(line[36:39])
        self.thetasolvals = [None]*self.nthetasolutes
        line = fileobj.readline() # Dummy line
        for i in range(self.nthetasolutes) :
          self.thetasolvals[i] = (fileobj.readline().strip().split()[2])
        
      line = fileobj.readline()
      if line.startswith("  -") or line.startswith("RESULTS FILE") : break
    # Sum of contributions for the different internal and interaction energies
    for attr in ["internal_energies","interaction_energies"] :
      if not hasattr(self,attr) : continue
      dict = getattr(self,attr)
      for label in dict :
        dict[label].append(EnergyResults())  
        for e in dict[label][:-1] :
          dict[label][-1] = dict[label][-1]+e
        dict[label][-1].type = "SUM"

    return line

class ResultsFile :
  """
  Class to store a collection of results

  Attributes
  ----------
  filename : string
    the name of the file
  snapshots : list of SnapshotResult
    all the results
  series : SnapshotResult
    all the results in a single object
  """
  def __init__(self) :
    self.filename = None
    self.snapshots = []
    self.series = None
  def read(self,filename,skip=0,readmax=None) :
    """
    Read results from disc
   
    The filename parameter can be either a single filename, in
    case it is assumed that it is a ProtoMS3.0 results file, i.e.
    it contains several snapshots
    if it is a list of strings, two behaviors are supported:
      1) if the length of the list is 2, and the first filename
    is a dictionary, a glob is used to list all filenames in that
    directory that starts with the second filename in the list
      2) if the length is not two or you give two filenames,
    each of them are processed individually and treated as individual
    snapshots

    Parameters
    ----------
    filename : string or list of strings
      the name of the file to read
    skip : int, optional
      number of snapshot at the beginning to skip
    readmax : int, optional
      maximum number of snapshots to read
    """
    if isinstance(filename,basestring) :
      f = open(filename,"r")
      line = f.readline()
      nsnap = 0
      while True :
        while line and not line.startswith("RESULTS FILE") : line = f.readline()
        if line :
          nsnap += 1
          if nsnap > skip :
            self.snapshots.append(SnapshotResults())
            line = self.snapshots[-1].parse(f)
          else :
            line = f.readline()
          if readmax is not None and len(self.snapshots) == readmax  :  break
        else :
          break
      f.close()  
    else :   # Assumes it is a list and read each one of them
      if readmax is None : readmax = 1E10
      if len(filename) == 2 :
        if os.path.isdir(filename[0]) :
          filenames = glob.glob(os.path.join(filename[0],"%s*"%filename[1]))
          if len(filenames) > 1 : filenames.sort()
        else :
          filenames = filename
      else :
        filenames = filename
      for filenam in filenames[skip:min(len(filenames),(skip+readmax+1))] :
        f = open(filenam, "r")
        self.snapshots.append(SnapshotResults())
        line = self.snapshots[-1].parse(f)
        f.close()
  def make_series(self) :
    """
    Make a snapshot with all the results

    This routine creates a special SnapshotResults object, which has
    the same attributes as the objects in the snapshots list,
    but rather than storing single floats and integers, this
    object stores NumpyArrays of all the results

    Returns
    -------
    SnapshotResults
      the created object, also stored in self.series
    """
    def set_energyresults(obj) :
      obj.curr = np.zeros(nsnap)
      obj.forw = np.zeros(nsnap)
      obj.back = np.zeros(nsnap)
    def put_energyresults(dest,source,i) :
      dest.curr[i] = source.curr
      dest.forw[i] = source.forw
      dest.back[i] = source.back
    
    if len(self.snapshots) < 2 : return
    nsnap = len(self.snapshots)

    # First replace all float/int in the SnapshotResults object with NumpyArrays
    self.series = copy.deepcopy(self.snapshots[0])
    for attr in ["lam","lamb","lamf","temperature","bvalues","pressure","volume","backfe","forwfe","gradient","agradient","capenergy"] :
      if hasattr(self.snapshots[0],attr) : setattr(self.series,attr,np.zeros(nsnap))
    for attr in ["datastep","lambdareplica","solventson","seed"] :
      if hasattr(self.snapshots[0],attr) : setattr(self.series,attr,np.zeros(nsnap,dtype=int))
    if hasattr(self.snapshots[0],"total") :
      set_energyresults(self.series.total)  
    if hasattr(self.snapshots[0],"internal_energies") :
      for elabel in self.series.internal_energies :
        for ene in self.series.internal_energies[elabel] :
          set_energyresults(ene)
    if hasattr(self.snapshots[0],"interaction_energies") :
      for elabel in self.series.interaction_energies :
        for ene in self.series.interaction_energies[elabel] :
          set_energyresults(ene)  
    if hasattr(self.snapshots[0],"extraenergy") :  
      set_energyresults(self.series.extraenergy)
    if hasattr(self.snapshots[0],"feenergies") :
      for elabel in self.series.feenergies :
        self.series.feenergies[elabel] = np.zeros(nsnap)
    if hasattr(self.snapshots[0],"thetavals") :
      for elabel in self.series.thetavals :
        self.series.thetavals[elabel] = np.zeros(nsnap)

    # Then loop over all snapshots and fill the NumpyArrays with data
    for i,snapshot in enumerate(self.snapshots) :
      for attr in ["lam","lamb","lamf","temperature","bvalues","pressure","volume","backfe","forwfe","gradient","agradient","datastep","lambdareplica","solventson","seed","capenergy"] :
        if hasattr(snapshot,attr) : getattr(self.series,attr)[i] = getattr(snapshot,attr)
      if hasattr(snapshot,"total") :
        put_energyresults(self.series.total,snapshot.total,i)
      if hasattr(snapshot,"internal_energies") :
        for elabel in snapshot.internal_energies :
          if elabel not in self.series.internal_energies : continue
          for ene1,ene2 in zip(self.series.internal_energies[elabel],snapshot.internal_energies[elabel]) :
             put_energyresults(ene1,ene2,i)
      if hasattr(snapshot,"interaction_energies") :
        for elabel in snapshot.interaction_energies :
          if elabel not in self.series.interaction_energies : continue
          for ene1,ene2 in zip(self.series.interaction_energies[elabel],snapshot.interaction_energies[elabel]) :
            put_energyresults(ene1,ene2,i)
      if hasattr(snapshot,"extraenergy") :  
        put_energyresults(self.series.extraenergy,snapshot.extraenergy,i)
      if hasattr(snapshot,"feenergies") :
        for elabel in snapshot.feenergies :
          if elabel not in self.series.feenergies : continue
          self.series.feenergies[elabel][i] = snapshot.feenergies[elabel]
      if hasattr(snapshot,"thetavals") :
        for elabel in snapshot.thetavals :
          if elabel not in self.series.thetavals : continue
          self.series.thetavals[elabel][i] = snapshot.thetavals[elabel]

    return self.series

#---------------------------------------------------------
# Classes to read, modify and write ProtoMS template files
#---------------------------------------------------------

class ForceFieldParameter :
  """ 
  Class to hold a general parameter set
  
  Attributes
  ----------
  key : string
    a label
  index : integer
    a serial number
  params : list
    either a list of strings or
    a list of ForceFieldParameter, containing
    the actual parameters
  """
  def __init__(self,record=None) :
    self.key = None
    self.index = None
    self.params = []
    if record is not None : self.parse(record)
  def parse(self,record) :
    """
    Parse a line from a ProtoMS file
    
    Parameters
    ----------
    record : string
      the line to parse from
    """
    cols = record.strip().split()
    self.key = cols[0] 
    self.index = int(cols[1])
    self.params = cols[2:]
  def __str__(self) :
    """
    Produces a correct ProtoMS file line
    """
    outparams = []
    for par in self.params :
      if isinstance(par,ForceFieldParameter) :
        outparams.append("%d"%par.index)
      else :
        outparams.append(par)
    return "%s %d    %s"%(self.key,self.index,"    ".join(outparams))

class AtomSet :
  """
  Class to hold a set of atoms and their associated parameter

  Attributes
  ----------
  atoms : list of strings
    the atoms in this set
  param : int or ForceFieldParameter
    the parameter associated with the atoms
  """
  def __init__(self,record=None) :
    self.atoms = []
    self.param = None
    if record is not None : self.parse(record)
  def parse(self,record) :
    """
    Parse a line from a ProtoMS file
    
    Parameters
    ----------
    record : string
      the line to parse from
    """
    cols = record.strip().split()
    self.atoms = cols[1:-1]
    self.param = int(cols[-1])
  def __str__(self) :
    """
    Produces a correct ProtoMS file line
    """
    if isinstance(self.param,ForceFieldParameter) :
      pindex = self.param.index
    else :
      pindex = self.param
    return "atm %s    %d"%(" ".join(self.atoms),pindex)

class TemplateAtom :
  """ 
  Class to hold a template atom

  Attributes
  ----------
  name : string
    the name of the atom
  residue : string
    the residue of the atom
  param0 : int or ForceFieldParameter
    the parameter associated with this atom at v0 
  param1 : int or ForceFieldParameter
    the parameter associated with this atom at v1
  bondedto : string or TemplateAtom
    the atom this atom is bonded to
  angleto : string or TemplateAtom
    the atom this atom forms an angle with
  dihedto : string or TemplateAtom
    the atom this atom forms a dihedral with
  """
  def __init__(self,record) :
    self.name = ""
    self.residue = ""
    self.param0 = None
    self.param1 = None
    self.bondedto = None
    self.angleto = None
    self.dihedto = None
    if record is not None : self.parse(record)
  def parse(self,record):
    """
    Parse a line from a ProtoMS file
    
    Parameters
    ----------
    record : string
      the line to parse from
    """
    pass
  def zmat(self) :
    """
    Produces a simple z-matrix representation of this atom
    
    Returns
    -------
    returns
      the z-matrix representation
    """
    pass

class TemplateSoluteAtom(TemplateAtom) :
  """ 
  Class to hold a solute atom

  This is a sub-class that implements functions
  that are specific for solute templates

  Has the same attributes as TemplateAtom

  """
  def __init__(self,record=None) :
    TemplateAtom.__init__(self,record)
  def parse(self,record) :
    """
    Parse a line from a ProtoMS file
    
    Parameters
    ----------
    record : string
      the line to parse from
    """
    cols = record.strip().split()
    self.name = cols[1] 
    self.residue = cols[2]
    self.param0 = int(cols[3])
    self.param1 = int(cols[4])
    self.bondedto = cols[5]
    self.angleto = cols[7]
    self.dihedto = cols[9]
  def __str__(self) :
    """
    Produces a correct ProtoMS file line
    """
    def make_str(strobj) :
      if isinstance(strobj,basestring) :
        if strobj in ["DM3","DM2","DM1"] :
          return "%4s %s"%(strobj,"DUM")
        else :
          return "%4s %s"%(strobj,self.residue)
      else :
        return "%4s %s"%(strobj.name,strobj.residue)
    def make_paramstr(param) :
      if isinstance(param,ForceFieldParameter) :
        return param.index
      else :
        return param
    return "atom %4s %s %d %d %s %s %s"%(self.name,self.residue,make_paramstr(self.param0),make_paramstr(self.param1),make_str(self.bondedto),make_str(self.angleto),make_str(self.dihedto))
  def zmat(self) :
    """
    Produces a simple z-matrix representation of this atom
    
    Returns
    -------
    returns
      the z-matrix representation, viz. ATOM BONDEDTO ANGLETO DIHEDRALTO
    """
    def make_str(strobj) :
      if isinstance(strobj,basestring) :
        return strobj
      else :
        return strobj.name  
    return "%s %s %s %s"%(self.name,make_str(self.bondedto),make_str(self.angleto),make_str(self.dihedto))
    
class TemplateConnectivity() :
  """ 
  Class to hold a solute bond, angle or dihedral

  Attributes
  ----------
  type : string
    the kind of connectivity, i.e. bond, angle or dihedral
  atoms : list of strings or TemplateAtom objects
    the atoms involved in this connectivity
  flex : float
    the flexibility
  param0 : int or ForceFieldParameter
    the parameter associated with this connectivity at v0 
  param1 : int or ForceFieldParameter
    the parameter associated with this connectivity at v1
  dummy : bool 
    flag to indicate if this should be treated as a dummy
  """
  def __init__(self,record=None) :
    self.type = ""
    self.atoms = []
    self.residues = []
    self.flex = None
    self.param0 = None
    self.param1 = None
    self.dummy = None
    if record is not None : self.parse(record)
  def parse(self,record) :
    """
    Parse a line from a ProtoMS file
    
    Parameters
    ----------
    record : string
      the line to parse from
    """
    cols = record.strip().split()
    self.type = cols[0]
    if self.type == "bond" :
      self.atoms = [cols[1],cols[3]]
      self.residues = [cols[2],cols[4]]
      nexti = 5
    elif self.type == "angle" :
      self.atoms = [cols[1],cols[3],cols[5]]
      self.residues = [cols[2],cols[4],cols[6]]
      nexti = 7
    elif self.type == "dihedral" :
      self.atoms = [cols[1],cols[3],cols[5],cols[7]]
      self.residues = [cols[2],cols[4],cols[6],cols[8]]
      nexti = 9
    while nexti < len(cols) :
      if cols[nexti] == "flex" :
        self.flex = float(cols[nexti+1])
        nexti = nexti + 2
      elif cols[nexti] == "param" :
        self.param0 = int(cols[nexti+1])
        self.param1 = int(cols[nexti+2])
        nexti = nexti + 3
      elif cols[nexti] == "dummy" :
        self.dummy = True
        nexti = nexti + 1
  def __str__(self) :
    """
    Produces a correct ProtoMS file line
    """
    def make_paramstr(param) :
      if isinstance(param,ForceFieldParameter) :
        return param.index
      else :
        return param
    def make_atomstr(atom) :
      if isinstance(atom,basestring) :
        return atom
      else :
        return atom.name
    strout = "%s %s"%(self.type," ".join("%4s %s"%(make_atomstr(atm),res) for atm,res in zip(self.atoms,self.residues)))
    if self.flex is not None : strout = strout + " flex %.3f"%self.flex          
    if self.param0 is not None : strout = strout + " param %d %d"%(make_paramstr(self.param0),make_paramstr(self.param1))
    if self.dummy : strout = strout + " dummy"
    return strout

class MolTemplate :
  """ 
  Class to hold a ProtoMS template  

  Attributes
  ----------
  name : string
    the name of the template
  type : string
    the kind of template, e.g. solute
  translate : float
    the translation displacement
  rot : float
    the rotational displacement  
  atoms : list of TemplateAtom objects
    the atoms in this template
  connectivity : list of TemplateConnectivity objects
    the connectivity in this template
  variables : list of string
    the variable geometries
  atomclass : TemplateAtom class
    the class to create objects of atoms
  """
  def __init__(self) :
    self.name = ""
    self.type = ""
    self.translate = 0.0
    self.rotate = 0.0
    self.jtheta = None
    self.jcorr = None
    self.jpmf = None
    self.atoms = []
    self.connectivity = []
    self.variables = []
    self.atomclass = None
  def parse_from(self,fileobj) :
    """
    Parse a full template from a file
    
    Terminates when found another mode 
    or end of file

    Parameters
    ----------
    fileobj : file object
      the file to parse from
    """
    line = fileobj.readline()
    while line :
      if line.startswith("solute") :
        self.type,self.name = line.strip().split()        
        self.atomclass = TemplateSoluteAtom
      elif line.startswith("info") :
        cols = line.strip().split()
        self.translate = float(cols[2])
        self.rotate = float(cols[4])
      elif line.startswith("jtheta") :
        self.jtheta = float(line.strip().split()[1])
      elif line.startswith("jcorr") :
        self.jcorr = float(line.strip().split()[1])
      elif line.startswith("jpmf") :
        self.jpmf = map(float,line.strip().split()[1:])
      elif line.startswith("atom") :
        self.atoms.append(self.atomclass(record=line))
      elif line.startswith("bond") or line.startswith("angle") or line.startswith("dihedral")  :
        self.connectivity.append(TemplateConnectivity(record=line))
      elif line.startswith("variable") :
        self.variables.append(line)
      elif line.startswith("mode") : return line
      line = fileobj.readline()
    return line
  def write_to(self,fileobj) :
    """
    Write the template to disc
   
    Parameters
    ----------
    fileobj : file object
      the file to write to
    """
    fileobj.write("mode template\n")
    fileobj.write("%s %s\n"%(self.type,self.name))
    fileobj.write("info translate %.3f rotate %.3f\n"%(self.translate,self.rotate))
    if self.jtheta is not None :
      fileobj.write("jtheta %.3f\n"%self.jtheta)
    if self.jcorr is not None :
      fileobj.write("jcorr %.3f\n"%self.jcorr)
    if self.jpmf is not None :
      fileobj.write("jpmf %s\n"%(" ".join("%.4f"%p for p in self.jpmf)))
    for atom in self.atoms : fileobj.write("%s\n"%atom)
    for con in self.connectivity : 
      # This was added to prevent angles and dihedral with dummy parameters from being sampled
      if len(con.atoms) > 2 and con.flex is not None and con.param0 == 0 and con.param1 == 0 : fileobj.write("#")
      fileobj.write("%s\n"%con)
    for var in self.variables : fileobj.write("%s\n"%var)
  def write_zmat(self,filename) :
    """
    Write the z-matrix of this template to disc
    
    Parameters
    ----------
    filename : string
      the name of the file to write to
    """
    with open(filename,"w") as f :
      for atom in self.atoms :
        f.write("%s\n"%atom.zmat())

class TemplateFile() :
  """ 
  Class to hold a ProtoMS template file

  Attributes
  ----------
  bondparams : list of ForceFieldParameter objects
    all bond parameters
  bondatoms : list of AtomSet objects
    all atoms associated with bond parameters
  angleparams : list of ForceFieldParameter objects
    all angle parameters
  angleatoms : list of AtomSet objects
    all atoms associated with angle parameters
  dihedralterms : list of ForceFieldParameter objects
    all dihedral term parameters
  dihedralparams : list of ForceFieldParameter objects
    all dihedral parameters
  dihedralatoms : list of AtomSet objects
    all atoms associated with dihedral parameters  
  cljparams : list of ForceFieldParameter objects
    all clj parameters
  templates : list of MolTemplate objects
    all templates in this file
  """   
  def __init__(self,filename=None) :
    self.bondparams = []
    self.bondatoms = []
    self.angleparams = []
    self.angleatoms = []
    self.dihedralterms = []
    self.dihedralparams = []
    self.dihedralatoms = []
    self.cljparams = []
    self.templates = []
    self.filename = ""
    if filename is not None : self.read(filename)

  def __str__(self) :
   return self.filename

  def append(self,other) :
    """
    Adds another template file to this one

    Will renumber the parameters of the other
    file
 
    Parameters
    ----------
    other : TemplateFile
      the file to add
    """
    templatenames = [t.name for t in self.templates]
    if set([t.name for t in other.templates]).issubset(set(templatenames)) :
      SetupError("All molecules in the new template file are already in the original template.")
      return

    if self.bondparams :
      start = max(p.index for p in self.bondparams) + 1
    else :
      start = -5000000
    for index,param in enumerate(other.bondparams,start) :
      if index > 0 : param.index = index
      self.bondparams.append(param)
    for atms in other.bondatoms :
      self.bondatoms.append(atms)

    if self.angleparams :
      start = max(p.index for p in self.angleparams) + 1
    else :
      start = -5000000
    for index,param in enumerate(other.angleparams,start) :
      if index > 0 : param.index = index
      self.angleparams.append(param)
    for atms in other.angleatoms :
      self.angleatoms.append(atms)

    if self.dihedralterms :
      start = max(p.index for p in self.dihedralterms) + 1
    else :
      start = -5000000
    for index,param in enumerate(other.dihedralterms,start) :
      if index > 0 : param.index = index
      self.dihedralterms.append(param)
    if self.dihedralparams :
      start = max(p.index for p in self.dihedralparams) + 1
    else :
      start = -5000000
    for index,param in enumerate(other.dihedralparams,start) :
      if index > 0 : param.index = index
      self.dihedralparams.append(param)
    for atms in other.dihedralatoms :
      self.dihedralatoms.append(atms)    

    if self.cljparams :
      start = max(p.index for p in self.cljparams) + 1    
    else :
      start = -50000
    for index,param in enumerate(other.cljparams,start) :
      if index > 0 : param.index = index
      self.cljparams.append(param)
    
    for template in other.templates :
      if template.name in templatenames :
        SetupError("Appending this template will cause duplicate names: %s. Aborting."%template.name)
      else :
        self.templates.append(template)     
     
  def assign_paramobj(self) :
    """
    Replaces integers and strings with objects

    It replaces all indices to force field parameters
    with references to ForceFieldParameter objects
 
    It replaces all atom in templates with references
    to TemplateAtom objects
    """
    def assign_con(con,paramlist,tem) :
      if not isinstance(con.param0,ForceFieldParameter) :
        for param in paramlist :
          if param.index == con.param0 :
            con.param0 = param
            break
      if not isinstance(con.param1,ForceFieldParameter) :
        for param in paramlist :
          if param.index == con.param1 :
            con.param1 = param
            break
      for i,atom in enumerate(con.atoms) :
        if isinstance(atom,basestring) :
          for atom2 in tem.atoms :
            if atom2.name.upper() == atom.upper() : 
              con.atoms[i] = atom2

    for atms in self.bondatoms :
      if not isinstance(atms.param,ForceFieldParameter) :
        for param in self.bondparams :
          if param.index == atms.param : 
            atms.param = param
            break
    for atms in self.angleatoms :
      if not isinstance(atms.param,ForceFieldParameter) :
        for param in self.angleparams :
          if param.index == atms.param : 
            atms.param = param
            break
    for param in self.dihedralparams :
      for i,term in enumerate(param.params) :
        if not isinstance(term,ForceFieldParameter) :
          for termff in self.dihedralterms : 
            if termff.index == int(term) :
              param.params[i] = termff
              break
    for atms in self.dihedralatoms :
      if not isinstance(atms.param,ForceFieldParameter) :
        for param in self.dihedralparams :
          if param.index == atms.param :
            atms.param = param
    for template in self.templates :
      for atom in template.atoms :
        if not isinstance(atom.param0,ForceFieldParameter) :
          for param in self.cljparams :
            if param.index == atom.param0 :
              atom.param0 = param
              break  
        if not isinstance(atom.param1,ForceFieldParameter) :
          for param in self.cljparams :
            if param.index == atom.param1 :
              atom.param1 = param
              break
        for atom2 in template.atoms :
          if isinstance(atom.bondedto,basestring) and atom.bondedto == atom2.name :
            atom.bondedto = atom2
          if isinstance(atom.angleto,basestring) and atom.angleto == atom2.name :
            atom.angleto = atom2
          if isinstance(atom.dihedto,basestring) and atom.dihedto == atom2.name :
            atom.dihedto = atom2
      for con in template.connectivity :
        if con.type == "bond" :
          assign_con(con,self.bondparams,template)
        elif con.type == "angle" :
          assign_con(con,self.angleparams,template)
        if con.type == "dihedral" :
          assign_con(con,self.dihedralparams,template)

  def read(self,filename) :
    """
    Read a template file from disc
 
    Parameters
    ----------
    filename : string
      the name of the file to read
    """    
    with open(filename,"r") as f :
      line = f.readline()
      while line :
        if line.startswith("mode bond") :
          line = f.readline()
          while line :
            if line.startswith("par") : 
              self.bondparams.append(ForceFieldParameter(record=line))      
            elif line.startswith("atm") : 
              self.bondatoms.append(AtomSet(record=line))
            elif line.startswith("mode") :
#              print "Found %s bondparams"%len(self.bondparams)
              break
            line = f.readline()
        elif line.startswith("mode angle") :
          line = f.readline()
          while line :
            if line.startswith("par") : 
              self.angleparams.append(ForceFieldParameter(record=line))      
            elif line.startswith("atm") : 
              self.angleatoms.append(AtomSet(record=line))
            elif line.startswith("mode") :
#              print "Found %s angleparams"%len(self.angleparams)
              break
            line = f.readline()
        elif line.startswith("mode dihedral") :
          line = f.readline()
          while line :
            if line.startswith("term") : 
              self.dihedralterms.append(ForceFieldParameter(record=line))      
            elif line.startswith("par") : 
              self.dihedralparams.append(ForceFieldParameter(record=line))
            elif line.startswith("atm") :
              self.dihedralatoms.append(AtomSet(record=line))
            elif line.startswith("mode") :
#              print "Found %s dihedralparams"%len(self.dihedralterms)
              break
            line = f.readline()
        elif line.startswith("mode clj") :
          line = f.readline()
          while line :
            if line.startswith("par") :
              self.cljparams.append(ForceFieldParameter(record=line))
            elif line.startswith("mode") :
              break
            line = f.readline()
        elif line.startswith("mode template") :
          self.templates.append(MolTemplate())
          line = self.templates[-1].parse_from(f)
    self.assign_paramobj()
    self.filename = filename

  def write(self,filename) :
    """
    Write a template file to disc
 
    Parameters
    ----------
    filename : string
      the name of the file to write to
    """    
    with open(filename,"w") as f :      
      if self.bondparams :
        f.write("mode bond\n")
        f.write("# U(r) = k(r-r0)**2\n")
        f.write("#parameter k(kcal mol-1 A-2) r0(A) comment\n")
        for param in self.bondparams : f.write("%s\n"%param)
        f.write("#atm atm1 atm2 parameter\n")
        for atms in self.bondatoms : f.write("%s\n"%atms)
      if self.angleparams  :
        f.write("mode angle\n")
        f.write("# U(theta) = k(theta-theta0)**2\n")
        f.write("#parameter k(kcal mol-1 deg-2) theta0(deg) comment\n")
        for param in self.angleparams : f.write("%s\n"%param)
        f.write("#atm atm1 atm2 atm3 parameter\n")
        for atms in self.angleatoms : f.write("%s\n"%atms)
      if self.dihedralparams  :
        f.write("mode dihedral\n")
        f.write("# U(phi) = k1*( 1.0 + k2*(cos[k3*phi + k4]) )\n")
        f.write("#term k1(kcal mol-1) k2 k3 k4(deg) #comment\n")
        for term in self.dihedralterms : f.write("%s\n"%term)
        f.write("#par  term1  term2 etc..  #comment\n")
        for par in self.dihedralparams : f.write("%s\n"%par)        
        f.write("#atm atm1 atm2 atm3 atm4 parameter #comment\n")
        for atms in self.dihedralatoms : f.write("%s\n"%atms)    
      if self.cljparams :
        f.write("mode clj\n")
        f.write("#parameter atm proton-num charge(|e|) sigma(A) epsilon(kcal mol-1) \n")
        for param in self.cljparams : f.write("%s\n"%param)
      if self.templates :
        for template in self.templates : template.write_to(f)

#----------------------------
# Classes to extend Argparse
#----------------------------

class _FullHelpAction(argparse.Action):
    """
    Class to initiate full help action
    """
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help=None):
        super(_FullHelpAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_fullhelp()
        parser.exit()

class MyArgumentParser(argparse.ArgumentParser) :
    """
    Extention to ArgumentParser to allow separation of of help in to 
    help for most common options and full help that list all arguments
    """
    def __init__(self,
                 prog=None,
                 usage=None,
                 description=None,
                 epilog=None,
                 version=None,
                 parents=[],
                 formatter_class=argparse.HelpFormatter,
                 prefix_chars='-',
                 fromfile_prefix_chars=None,
                 argument_default=None,
                 conflict_handler='error',
                 add_help=True,
                 add_fullhelp=True):
        super(MyArgumentParser,self).__init__(
                  prog,usage,description,epilog,version,parents,formatter_class,
                  prefix_chars,fromfile_prefix_chars,argument_default,
                  conflict_handler,add_help)
               
        self.add_fullhelp = add_fullhelp          
        self.register('action', 'fullhelp', _FullHelpAction)
                  
        if '-' in prefix_chars:
            default_prefix = '-'
        else:
            default_prefix = prefix_chars[0]
        if self.add_fullhelp:
            self.add_argument(
                default_prefix*2+'fullhelp',
                action='fullhelp', default=argparse.SUPPRESS,
                help=gettext.gettext('show full help and exit'))

        self._action_groups[1].title = "Most common arguments"

    def format_help(self,fullhelp=False):
        formatter = self._get_formatter()

        # usage
        formatter.add_usage(self.usage, self._actions,
                            self._mutually_exclusive_groups)

        # description
        formatter.add_text(self.description)

        # positionals, optionals and user-defined groups
        if fullhelp :
          grps = self._action_groups
        else :
          grps = self._action_groups[:2]
        for action_group in grps :
            formatter.start_section(action_group.title)
            formatter.add_text(action_group.description)
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        # epilog
        formatter.add_text(self.epilog)
        
        if not fullhelp :  
          formatter.add_text("Type %s --fullhelp for complete list of all arguments"%self.prog)

        # determine help from format above
        return formatter.format_help()

    def print_help(self, file=None):
        if file is None:
            file = sys.stdout
        self._print_message(self.format_help(fullhelp=False), file)

    def print_fullhelp(self, file=None):
        if file is None:
            file = sys.stdout
        self._print_message(self.format_help(fullhelp=True), file)

#-----------------------
# Other useful routines
#-----------------------

def standard_filename(filename,folder) :
    """ 
    Generates the filename of file in the standard ProtoMS file hierarchy

    If $PROTOMSHOME is set, it uses it as the base-path, otherwise,
    it uses the location of this module as to create the base-path

    Parameters
    ----------
    filename : string
      the name of the file
    folder : string
      a folder in the standard ProtoMS file hierarchy
   
    Returns
    -------
    the full path to the file
    """
    folder_filename = os.path.join(folder,filename)
    pmshome = os.getenv("PROTOMSHOME")
    # If $PROTOMSHOME is set, use that as the base path
    if pmshome is not None :
      return os.path.join(pmshome,folder_filename)
    else :
    # otherwise, use the location of this python library,
    # and step-up onces in the hierarchy to find the base path
      thispath = os.path.dirname(os.path.abspath(__file__))
      oneup = os.path.split(thispath)[0]
      return os.path.join(oneup,folder_filename)

def setup_logger(filename=None) :
  """
  Setup ProtoMS logging system
  
  Uses the standard logging module with two handlers,
  one that prints everything above DEBUG to standard output
  and one that prints everything, even DEBUG to a file
  
  If filename is None no file handler is created
  
  This should be called by protoms.py and all other stand-alone
  programs that uses the tools
  
  Parameters
  ----------
  filename : string, optional
    the filename of the string
    
  Returns
  -------
  a reference to the created logger
  """
  logger = logging.getLogger('protoms')
  logger.setLevel(logging.DEBUG)
  formatter1 = logging.Formatter('%(message)s')
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  console.setFormatter(formatter1)
  logger.addHandler(console)
  if filename is not None :
    formatter2 = logging.Formatter('%(levelname)s : %(message)s')
    logfile = logging.FileHandler(filename,mode="w")
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(formatter2)
    logger.addHandler(logfile)
  return logger
  
def angle(v1,v2) :
  """
  Calculates the angle between two vectors
  
  Parameters
  ----------
  v1 : Numpy array
    the first vecor
  v2 : Numpy array
    the second vecor

  Returns
  -------
  float
    the angle in radians
  """
  a = (v1*v2).sum()/(np.sqrt((v1**2).sum())*np.sqrt((v2**2).sum()))
  if a > 1.0 :
    a = 0.0
  elif a < -1 :
    a = np.pi
  else :
    a = np.arccos(a)
  return a

def angle_atms(a1,a2,a3) :
  """
  Calculates the angle between three atoms
  
  Parameters
  ----------
  a1 : Numpy array
    the first atom
  a2 : Numpy array
    the second atom
  a3 : Numpy array
    the third atom

  Returns
  -------
  float
    the angle in radians
  """
  v1 = a2 - a1
  v2 = a2 - a3
  return angle(v1,v2)

def color(idx) :
  """
  Returns a color of index

  For instances when the index is larger than the number of defined colors,
  this routine takes care of this by periodicity, i.e.
  color at idx=0 is the same color as idx=n 

  Parameters
  ----------
  idx : int
    the index of the color
  
  Returns
  -------
  list of floats 
    the color
  """
  colors = []
  colors.append((0.0/255.0,69.0/255.0,134.0/255.0))
  colors.append((255.0/255.0,66.0/255.0,14.0/255.0))
  colors.append((255.0/255.0,211.0/255.0,32.0/255.0))
  colors.append((87.0/255.0,157.0/255.0,28.0/255.0))
  colors.append((126.0/255.0,0.0/255.0,33.0/255.0))
  colors.append((131.0/255.0,202.0/255.0,255.0/255.0))
  colors.append((49.0/255.0,64.0/255.0,4.0/255.0))
  colors.append((174.0/255.0,207.0/255.0,0.0/255.0))
  d = int(len(colors)*np.floor(idx / float(len(colors))))
  return colors[idx-d]

if __name__ == "__main__":

  raise SetupError("This module cannot be executed as a stand-alone program!")
        
