# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Julien Michel
#          Gregory Ross

""" Useful classes for setup and analysis of ProtoMS simulations.
    Contain classes to 
      1) Read, modify and write structures in PDB format
      2) Store parameter collections
      3) Read and store ProtoMS results files
      4) Read, modify and write ProtoMS template files
"""

import math
import sys
import os
import numpy as np

boltz = 0.00198717076 # kcal.mol-1K-1

class SetupError(Exception) :
    """ A general exception class to be raised by setup code
    """
    pass

def standard_filename(filename,folder) :
    """ Returns the filename of file in the
        standard ProtoMS file hierarchy
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

#------------------------------------------
# Classes and routines to handle PDB-files
#------------------------------------------

class Atom(object):
    """ Class for holding a PDB atom record
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
        """ Trying to set the element of this atom
        """
        k = 0
        self.element = self.name[k]
        while self.element.isdigit():
            k = k + 1
            self.element = self.name[k]
        self.element = self.element.lower()
    def __str__(self):
        return "ATOM  %-5d %3s %3s    %4d     %8.5f %8.5f %8.5f  1.00 0.00" % (self.index,self.name,self.resname,self.resindex,self.coords[0],self.coords[1],self.coords[2])
              
class Residue(object):
    """ Class for holding a set of PDB atoms
    """
    def __init__(self,name="???",index=0):
        self.atoms = []
        self.index = index
        self.name = name
        self.center = np.array([0.0,0.0,0.0])          
    def addAtom(self,atom=None):
        """ Adds an atom to this residue
        """        
        assert atom.type ==  "atom"
        self.atoms.append(atom)
    def getCenter ( self ):
        """ Calculates the center of coordinates for this residue 
        """
        coords = np.array ( [ atom.coords for atom in self.atoms ] )
        self.center = coords[:,0].mean(), coords[:,1].mean(), coords[:,2].mean()
        #return self.center            
    def __str__(self):
        return '\n'.join(atom.__str__() for atom in self.atoms)

class PDBFile:
    """ Class for holding the atoms and residues of a PDB file
    """
    def __init__ ( self, filename = None ):
        self.residues = {}
        self.solvents = {}
        self.name = ""
        self.center = np.array([0.0,0.0,0.0])
        self.header = ""

        if filename is not None : self.read(filename)

    def __str__ ( self ):
        return self.name

    def copy ( self ) :
        """ Make a copy of the residues and solvents dictionaries, not the Residue and Atom objects themselves
        """
        new = PDBFile()
        new.name = self.name
        new.header = self.header
        new.center = np.array(self.center,copy=True)
        for res in self.residues : new.residues[res] = self.residues[res]
        for sol in self.solvents : new.solvents[sol] = self.solvents[sol]
        return new

    def read(self, filename) :
        residues = {}
        solvents = {}
        self.name = filename
        with open ( filename ) as f:
          line = f.readline()
          nres = 0
          prevres = -1
          while line :
            if line[:6] in ["ATOM  ","HETATM"] :                
                index = int(line[6:11].strip())
                atname = line[12:16].strip()
                restype = line[17:20]
                resnum = int(line[22:26].strip())
                if resnum != prevres :
                  nres = nres + 1
                prevres = resnum
                x,y,z = float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())
                coords = [x,y,z]
                newatom = Atom(index=index,name=atname,resindex=nres,
                               resname=restype,coords=coords)
                # If solvent
                if restype not in ['WAT','HOH','T3P','t3p','T4P','t4p','SOL','sol','see','SEE']:
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
            # Read next line
            line = f.readline()
        self.residues, self.solvents = residues, solvents

    def write ( self, filename, renumber = False, header = None ):
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
            f.write ( "TER \n" )
            for i, sol in enumerate ( sorted ( self.solvents.keys() ), 1 ):
                for atom in self.solvents[sol].atoms:
                    if renumber:
                        atom.resindex = i
                    s = "ATOM  %5d %-4s %3s  %4d    %8.3f%8.3f%8.3f        \n" % (atom.index+1,atom.name,atom.resname,atom.resindex,atom.coords[0],atom.coords[1],atom.coords[2])
                    f.write ( s )
                f.write ( "TER \n" )

    def getCenter(self):
        center = np.array([0.0,0.0,0.0])
        for res in self.residues:
            self.residues[res].getCenter()
            rescent = self.residues[res].center
            center = center + rescent
        self.center = center/float(len(self.residues))
        

def merge_pdbs(pdbobjs) :
    pdbout = PDBFile()
    nres = 0
    for pdbobj in pdbobjs :
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
    return pdbout  
      
#--------------------------------
# Classes to hold parameter sets
#--------------------------------

class parameter:
    #Not valid for dihedral parameters but contains sufficient
    #information to check that a parameter exists
    def __init__ ( self, index, ats, k, b0 ):
        self.index = int ( index )
        self.ats = ats
        self.k = float ( k )
        self.b0 = float ( b0 )

class parameter_set:
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
            self.params += [ parameter ( par_cols[1], at_cols[1:-1], par_cols[2], par_cols[3] ) ]

    def get_params ( self, ats ):
        try:
            return [ i for i in self.params 
                     if i.ats == ats or i.ats == ats[::-1] ][0]
        except IndexError:
            return [ i for i in self.params 
                     if i.ats[1:3] in [ ats[1:3], ats[1:3:-1] ] ][0]

#--------------------------------------------------
# Classes to read and store ProtoMS results files
#--------------------------------------------------

class EnergyResults :
  def __init__(self,line=None) : 
    self.curr = None
    self.back = None
    self.forw = None
    self.type = None
    if line is not None :
      self.parse_line(line)
  def parse_line(self,line) :
    cols = line.strip().split()
    self.type = cols[0] 
    self.curr = float(cols[1])
    self.forw = float(cols[6])
    self.back = float(cols[10])
  def __str__(self) :
    return "%s %F20.10 %F20.10 %F20.10"%(self.type,self.curr,self.back,self.forw)

class SnapshotResults :
  def __init__(self,fileobj=None) :
    if fileobj is not None :
      self.parse(fileobj) 
  def parse(self,fileobj) :
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
      if line.startswith(" Simulation B factor") : self.bfactor = float(line.split("=")[1].strip())
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
    self.average_energies = {}
    while line :
      if line.startswith("Internal ")  :
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
          key = "protein"+cols[4]+"-solute"+cols[8]
        elif cols[1] == "protein-solvent" :
          key = "protein"+cols[4]+"-solvent"
        elif cols[1] == "solute-solute" :
          key = "solute"+cols[5]+"-solute"+cols[8]
        elif cols[1] == "solute-solvent" :
          key = "solute"+cols[5]+"-solvent"
        elif cols[1] == "solute-GCS" :
          key = "solute"+cols[5]+"-GCS"   
        self.average_energies[key] = []   
        line = fileobj.readline() # Dummy line
        self.average_energies[key].append(EnergyResults(line=fileobj.readline())) # Coul
        line = fileobj.readline() # Dummy line    
        self.average_energies[key].append(EnergyResults(line=fileobj.readline())) # LJ
      elif line.startswith("FREE ENERGY DATA") :
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # Lambda header
        line = fileobj.readline() # Dummy line
        line = fileobj.readline() # First lambda
        self.feenergies = {}
        while line[0] != "#" :
          cols = line.strip().split()
          self.feenergies[float(cols[0])] = float(cols[2])
          line = fileobj.readline()
        self.gradient = float((fileobj.readline().strip().split()[1]))
      elif line.startswith(" Printing individual") :
        self.thetavals = [None]*self.ngcsolutes
        line = fileobj.readline() # Dummy line
        for i in range(self.ngcsolutes) :
          self.thetavals[i] = (fileobj.readline().strip().split()[2])
        
      line = fileobj.readline()
      if line.startswith("  -") : break

class ResultsFile :
  def __init__(self) :
    self.filename = None
    self.snapshots = []
  def read(self,filename,skip=0,readmax=None) :
    if isinstance(filename,basestring) :
      f = open(filename,"r")
      line = f.readline()
      nsnap = 0
      while True :
        while line and not line.startswith("RESULTS FILE") : line = f.readline()
        if line :
          snapid = int(line.strip().split()[2])
          if snapid > skip :
            self.snapshots.append(SnapshotResults(f))
          if readmax is not None and len(self.snapshots) == readmax  : break
        else :
          break
        line = f.readline()
      f.close()  
    else :   # Assumes it is a list and read each one of them
      if readmax is None : readmax = 1E10
      for filenam in filename[skip:min(len(filename),(skip+readmax+1))] :
        f = open(filenam, "r")
        self.snapshots.append(SnapshotResults(f))
        f.close()
        
#---------------------------------------------------------
# Classes to read, modify and write ProtoMS template files
#---------------------------------------------------------

class ForceFieldParameter :
  """ Class to hold a general parameter set
  """
  def __init__(self,record=None) :
    self.key = None
    self.index = None
    self.params = []
    if record is not None : self.parse(record)
  def parse(self,record) :
    cols = record.strip().split()
    self.key = cols[0] 
    self.index = int(cols[1])
    self.params = cols[2:]
  def __str__(self) :
    outparams = []
    for par in self.params :
      if isinstance(par,ForceFieldParameter) :
        outparams.append("%d"%par.index)
      else :
        outparams.append(par)
    return "%s %d    %s"%(self.key,self.index,"    ".join(outparams))

class AtomSet :
  """ Class to hold a set of atoms and their associated parameter
  """
  def __init__(self,record=None) :
    self.atoms = []
    self.param = None
    if record is not None : self.parse(record)
  def parse(self,record) :
    cols = record.strip().split()
    self.atoms = cols[1:-1]
    self.param = int(cols[-1])
  def __str__(self) :
    if isinstance(self.param,ForceFieldParameter) :
      pindex = self.param.index
    else :
      pindex = self.param
    return "atm %s    %d"%(" ".join(self.atoms),pindex)

class TemplateAtom :
  """ Class to hold a template atom
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
    pass
  def zmat(self,record) :
    pass

class TemplateSoluteAtom(TemplateAtom) :
  """ Class to hold a solute atom
  """
  def __init__(self,record=None) :
    TemplateAtom.__init__(self,record)
  def parse(self,record) :
    cols = record.strip().split()
    self.name = cols[1] 
    self.residue = cols[2]
    self.param0 = int(cols[3])
    self.param1 = int(cols[4])
    self.bondedto = cols[5]
    self.angleto = cols[7]
    self.dihedto = cols[9]
  def __str__(self) :
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
    def make_str(strobj) :
      if isinstance(strobj,basestring) :
        return strobj
      else :
        return strobj.name  
    return "%s %s %s %s"%(self.name,make_str(self.bondedto),make_str(self.angleto),make_str(self.dihedto))
    
class TemplateConnectivity() :
  """ Class to hold a solute bond, angle or dihedral
  """
  def __init__(self,record=None) :
    self.type = ""
    self.atoms = []
    self.flex = None
    self.param0 = None
    self.param1 = None
    self.dummy = None
    if record is not None : self.parse(record)
  def parse(self,record) :
    cols = record.strip().split()
    self.type = cols[0]
    if self.type == "bond" :
      self.atoms = ["%4s %s"%(cols[1],cols[2]),"%4s %s"%(cols[3],cols[4])]
      nexti = 5
    elif self.type == "angle" :
      self.atoms = ["%4s %s"%(cols[1],cols[2]),"%4s %s"%(cols[3],cols[4]),"%4s %s"%(cols[5],cols[6])]
      nexti = 7
    elif self.type == "dihedral" :
      self.atoms = ["%4s %s"%(cols[1],cols[2]),"%4s %s"%(cols[3],cols[4]),"%4s %s"%(cols[5],cols[6]),"%4s %s"%(cols[7],cols[8])]
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
    def make_paramstr(param) :
      if isinstance(param,ForceFieldParameter) :
        return param.index
      else :
        return param
    strout = "%s %s"%(self.type," ".join(self.atoms))
    if self.flex is not None : strout = strout + " flex %.3f"%self.flex          
    if self.param0 is not None : strout = strout + " param %d %d"%(make_paramstr(self.param0),make_paramstr(self.param1))
    if self.dummy : strout = strout + " dummy"
    return strout

class MolTemplate :
  """ Class to hold a ProtoMS template  
  """
  def __init__(self) :
    self.name = ""
    self.type = ""
    self.translate = 0.0
    self.rotate = 0.0
    self.atoms = []
    self.connectivity = []
    self.variables = []
    self.atomclass = None
  def parse_from(self,fileobj) :
    line = fileobj.readline()
    while line :
      if line.startswith("solute") :
        self.type,self.name = line.strip().split()        
        self.atomclass = TemplateSoluteAtom
      elif line.startswith("info") :
        cols = line.strip().split()
        self.translate = float(cols[2])
        self.rotate = float(cols[4])
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
    fileobj.write("mode template\n")
    fileobj.write("%s %s\n"%(self.type,self.name))
    fileobj.write("info translate %.3f rotate %.3f\n"%(self.translate,self.rotate))
    for atom in self.atoms : fileobj.write("%s\n"%atom)
    for con in self.connectivity : fileobj.write("%s\n"%con)
    for var in self.variables : fileobj.write("%s\n"%var)
  def write_zmat(self,filename) :
    with open(filename,"w") as f :
      for atom in self.atoms :
        f.write("%s\n"%atom.zmat())

class TemplateFile() :
  """ Class to hold a ProtoMS template file
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
    if filename is not None : self.read(filename)

  def append(self,other) :

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
    
    templatenames = [t.name for t in self.templates]
    for template in other.templates :
      if template.name in templatenames :
        SimulationObjects.SetupError("Appending this template will cause duplicate names: %s. Aborting."%template.name)
      self.templates.append(template)     
     
  def assign_paramobj(self) :
    def assign_con(con,paramlist) :
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
      for con in template.connectivity :
        if con.type == "bond" and con.param0 is not None :
          assign_con(con,self.bondparams)
        elif con.type == "angle" and con.param0 is not None :
          assign_con(con,self.angleparams)
        if con.type == "dihedral" and con.param0 is not None :
          assign_con(con,self.dihedralparams)
  def read(self,filename) :
    
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

  def write(self,filename) :

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

if __name__ == "__main__":

  raise SetupError("This module cannot be executed as a stand-alone program!")
        
