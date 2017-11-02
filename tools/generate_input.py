# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Classes and routines to make ProtoMS command files

This module defines a single public function
generate_input

and the following public classes:
ProtoMSSimulation
ProteinLigandSimulation
Equilibration
Sampling
DualTopology

Can be executed from the command line as a stand-alone program
"""

import os, sys
import logging

import numpy as np

import random

import simulationobjects as sim
import make_single

logger = logging.getLogger('protoms')

# Default force constant for absolute binding free energies
FORCE_CONSTANT = 1

def _assignMoveProbabilities(protein,solute,solvent,moveset,isperiodic) :
  """ 
  Assigns move probabilities for protein, solute, solvent and box
  Does this by some "heuristic" rule
  
  Parameters
  ----------
  protein : string
    the filename of the protein pdb file
  solute : list of strings
    the filenames of the solute pdb file
  solvent : string
    the filename of the solvent pdb file
  moveset : string
    string indicating the type of move set to create.
    Can be standard, gcmc, jaws1 or jaws2
  isperiodic : boolean
    flag indicating if a periodic simulation is prepared
  """
  pvolume = 0.0
  numprot = 0.0
  numsolu = 0.0
  numsolv = 0.0
  pinsert = 0.0
  pdelete = 0.0
  pgcsolu = 0.0

  accel_prot = 5					# The amount by which the protein be accelerated over the solvent.
  accel_solu = 50       				# The amount by which the solute be accelerated over the solvent.
  addto = 1000.0					# Move proportion add up to this number. This is just for clarity.

  if protein is not None:
    protobj = sim.PDBFile(filename=protein)
    numprot = float(len(protobj.residues))*accel_prot
  if solute is not None:
    soluobj =  sim.PDBFile(filename=solute[0])
    numsolu = float(len(soluobj.residues))*accel_solu
  if solvent is not None:
    solvobj =  sim.PDBFile(filename=solvent)
    numsolv = float(len(solvobj.solvents))
  
  numtot = numprot + numsolu + numsolv
  if isperiodic:
    pvolume = int(0.002*addto) + 1
  else:
    pvolume = 0
  if moveset == "standard" :
    psolvent = int(round(addto*numsolv/numtot)) - pvolume
    psolute = int(round(addto*numsolu/numtot))
    pprotein = int(round(addto*numprot/numtot))
    return "solvent=%d protein=%d solute=%d volume=%d"%(psolvent,pprotein,psolute,pvolume)
  elif moveset == "gcmc":
    pinsert = int(round(addto/6))				# Giving half of all moves to GCMC moves.
    pdelete = int(round(addto/6))
    pgcsolu = int(round(addto/6))
    psolvent = int(round(addto*numsolv/numtot/2.0))
    psolute = int(round(addto*numsolu/numtot/2.0))
    pprotein = int(round(addto*numprot/numtot/2.0))
    return "solvent=%d protein=%d solute=%d insertion=%d deletion=%d gcsolute=%d"%(psolvent,pprotein,psolute,pinsert,pdelete,pgcsolu)
  elif moveset in ["jaws1","jaws2"] :
    pgcsolu = int(round(addto/6)) # Giving half of all move to JAWS moves.
    ptheta  = int(round(addto/3))
    psolvent = int(round(addto*numsolv/numtot/2.0))
    psolute = int(round(addto*numsolu/numtot/2.0))
    pprotein = int(round(addto*numprot/numtot/2.0))
    if moveset == "jaws1" :
      movenam = "theta"
    else :
      movenam = "sample"
    return "solvent=%d protein=%d solute=%d %s=%d gcsolute=%d"%(psolvent,pprotein,psolute,movenam,ptheta,pgcsolu)
  
class ProtoMSSimulation :
    """
    This is a ProtoMS command file object.This object holds all of the information
    about the simulation and can either run the simulation with this
    information, or can write this information to a command file.
       
    Attributes
    ----------
    lines : 
        list of strings all the lines of the ProtoMS command file    
    """
    def __init__(self,filename=None):
        """
         Parameters
        ----------
        filename : string, optional
            the filename of the command file to read from disc
        """
        self._lines = []
        if (not filename is None):
            self.readCommandFile(filename)        
    def __str__(self):
        string = "\n".join(self._lines)           
        return string        
    def _setkey(self,key,value):
      #print [key in line.split()[0] for line in self._lines]
      self._lines.append("%s %s"%(key,value))        
    def clear(self):
        """
        Clear the simulation description
        """    
        self._lines = []
    def readCommandFile(self,filename):
        """
        Read in a simulation from a command file
        
        Parameters
        ----------
        filename : string
            the filename of the command file to read from disc
        """
        lines = open(filename,"r").readlines()  
 
        for line in lines:
            #split the line into words
            words = line.split()
            #skip short lines or comments
            if (len(words) < 2): continue
            if (words[0].find("#") == 0): continue

            arguments = ' '.join(words[1:])
            self._setkey(words[0],arguments)
    def setStream(self,stream,filename):
        """
        Set the direction of the output stream to a filename
        
        Parameters
        ----------
        stream : string
          the name of the stream
        filename : string
          the filename to direct the stream to
        """
        key = "stream%s" % stream.lower()
        self._setkey(key,filename)
    def setParameter(self,parameter,value):
        """
        Set the parameter to value
        
        Parameters
        ----------
        parameter : string
          the name of the parameter
        value : string
          the value of the parameter
        """
        self._setkey(parameter.lower(),value)
    def setChunk(self,chunk):
        """
        This adds the simulation chunk to 'chunk'.
        
        Parameters
        ----------
        chunk : string
          the chunk to be executed, including parameters
        """
        self._setkey("chunk",chunk)
    def setDump(self,dump,freq):
        """
        This adds the simulation dump to 'dump' with frequnecy freq
        
        Parameters
        ----------
        dump : string
          the dump to be executed, including parameters
        freq : int
          the dump frequency
        """
        self._setkey("dump","%d %s"%(freq,dump))
    def setProtein(self,n,filename):
        """
        Set the filename for proteinN. It then returns n+1
        
        Parameters
        ----------
        n : int
          the protein serial number
        filename : string
          the filename of the protein pdb file
          
        Returns
        -------
        int
            n + 1
        """
        key = "protein%d" % n
        self._setkey(key,filename)
        return n+1      
    def setSolute(self,n,filename):
        """
        Set the filename for soluteN. It then returns n+1
                
        Parameters
        ----------
        n : int
          the solute serial number
        filename : string
          the filename of the solute pdb file
        
        Returns
        -------
        int
            n + 1
        """
        key = "solute%d" % n
        self._setkey(key,filename)
        return n+1
    def setSolvent(self,n,filename):
        """
        Set the filename for solventN. It then returns n+1
                
        Parameters
        ----------
        n : int
          the solvent serial number
        filename : string
          the filename of the solvent pdb file
        
        Returns
        -------
        int
            n + 1
        """
        key = "solvent%d" % n
        self._setkey(key,filename)
        return n+1
    def setForceField(self,filename):
        """
        Set the filename for parfile
        
        Parameters
        ----------
        filename : string
          the filename of the force field file
        """
        self._setkey("parfile",filename)
    def writeCommandFile(self,filename):
        """
        Write the contents of this simulation object to disc
        
        Parameters
        ----------
        filename : string
          the filename to write the commands to
        """
        f = open(filename,"w")
        for line in self._lines:
            f.write("%s\n" % line)            
        f.close()
    def run(self,exe,cmdfile=None):
        """
        Run this simulation
        
        Parameters
        ----------
        exe : string
          the filename of the executable
        cmdfile : string, optional
          the filename of the command file
        """
        if (cmdfile is None):
            #save commands as environmental variables
            for line in self._lines:
                words = line.split()
                os.putenv(words[0],join(word[1:]))         
            #run the executable
            exitval = subprocess.call(exe,shell=True)
        else:
            # save commands as a file
            writeCommandFile(cmdfile)
            exitval = subprocess.call("%s %s" % (exe,cmdfile),shell=True)
        return exitval

class ProteinLigandSimulation(ProtoMSSimulation) :
  """ 
  This a generic command file for a protein-ligand simulation
  Generates input for proteins, solutes and solvent but does not write any chunks
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    outfolder=None,
                    ranseed=None) :
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    outfolder : string
      the output folder for PROTOMS
    """
    
    ProtoMSSimulation.__init__(self)
    self.setForceField(os.path.join("$PROTOMSHOME","parameter","amber14SB.ff"))
    self.setForceField(os.path.join("$PROTOMSHOME","parameter","solvents.ff"))
    self.setForceField(os.path.join("$PROTOMSHOME","parameter","amber14SB-residues.ff"))
    self.setForceField(os.path.join("$PROTOMSHOME","parameter",ProteinLigandSimulation.gaffversion + ".ff"))
    if templates is not None and templates :
      for tem in templates :
        self.setForceField(tem)
    if not protein is None :
      self.setProtein(1,protein)
    if solutes :
      for i,sol in enumerate(solutes) :
        self.setSolute(i+1,sol)
    self.periodic = False
    if solvent is not None :
      self.setSolvent(1,solvent)
      try :
        for line in open(solvent) :
          if line.lower().find("header box") > -1 :
            self.periodic = True
            break 
      except :
        pass      
    if outfolder is not None: self.setParameter("outfolder",outfolder)
    self.setStream("header","off")
    self.setStream("detail","off")
    for s in ["warning","info","fatal","results","accept",] :
      self.setStream(s,s)
    self.setParameter("cutoff",10.0)
    self.setParameter("feather",0.5)
    self.setParameter("temperature",25.0)
    if ranseed is None:				
      r = random.randrange(110000,11000000)
      if r % 2 == 0 : r = r + -1
      self.setParameter("ranseed",r)
    else:
      self.setParameter("ranseed",ranseed)
    if not solvent is None :
      self.setParameter("boundary","solvent")
      if self.periodic  :
        self.setParameter("pressure","1")
    if protein is None :
      for i,sol in enumerate(solutes) :
        if solvent is None :
          self.setChunk("transrot %d 0.0 0.0"%(i+1))
        else :
          self.setChunk("transrot %d 0.0"%(i+1))
    else :
      self.setParameter("pdbparams","on")

class Equilibration(ProteinLigandSimulation) :
  """ 
  This a command file for equilibration of a protein-ligand system
  Generates input for proteins, solutes and solvent and write an equilibration and a pdb chunk
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    outfolder="",
                    ranseed=None,
                    nsteps=1E5,
                    pdbfile="equilibrated.pdb") :
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    outfolder : string
      the output folder for PROTOMS
    nsteps : int, optional
      number of equilibration moves
    pdbfile  : string, optional
      name of the filename to save the equilibrated system to
    """                
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)

    #moves = assignMoveProbabilities(protein is not None,solutes,solvent is not None,self.periodic)
    moves = _assignMoveProbabilities(protein,solutes,solvent,"standard",self.periodic)
    self.setChunk("equilibrate %d %s"%(nsteps,moves))
    self.setChunk("pdb all solvent=all file=%s standard"%pdbfile)

class Sampling(ProteinLigandSimulation) :
  """ 
  This a command file for MC sampling of a protein-ligand system
  Generates input for proteins, solutes and solvent and write 
  chunks for equilibration and simulation, as well as dumps
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    outfolder="",
                    ranseed=None,
                    nequil=5E6,
                    nprod=40E6,
                    dumpfreq=1E5,
                    outprefix="",
                    tune=False) :
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    outfolder : string
      the output folder for PROTOMS
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    outprefix  : string, optional
      the prefix for all output files
    tune : boolean, optional
      perform a dihedral tuning simulation, default False
    """                 
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)

    if tune:
      self.setParameter ( 'tunedihedral', 'on' )

    self.setDump("results write results",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"standard",self.periodic)
    if nequil > 0 : self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class DualTopology(ProteinLigandSimulation) :
  """ 
  This a command file for a dual-topology protein-ligand simulation
  Generates input for proteins, solutes and solvent
  write dual-toplogy parameters, lambda-replica exchange,
  dump results, pdb and restart files
  equilibrate and production run
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute1.pdb","solute2.pdb"],
                    solvent="water.pdb",
                    templates=["solute1.tem","solute2.tem"],
                    nequil=5E6,
                    nprod=40E6,
                    dumpfreq=1E5,
                    ranseed=None,
                    lambdaval=None,
                    outfolder="out",
                    restrained=[],
                    softcore='mixed') :  
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    lambdaval : float or list of floats
      the lambda values to perform the simulation at
    outfolder  : string, optional
      the folder for all output files
    restrained : string, optional
      the solutes on which restrains should be applied
    sofcore : string, optional
      determine which atoms to apply softcore potentials to
      can have value of 'all', 'none', 'mixed'
      'all' - softcores applied to all atoms of both solutes
      'none' - softcores not applied to any atoms
      'mixed' - softcores will be applied only to non matching 
                atoms within ligand structures
    
    Raises
    ------
    SetupError
      less than two solutes given or no lambda values given
    """      
    if len(outfolder) == 0 : outfolder = "out"          
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)

    if len(solutes) < 2 and restrained is [] :
      raise sim.SetupError("Cannot do dual topology with less than 2 solutes")
      
    if lambdaval is None or len(lambdaval) < 2 :
      raise sim.SetupError("Must give at least two lambda values")

    if softcore not in ('all', 'none', 'mixed'):
      raise sim.SetupError("softcore argument must be one of 'all', 'none' or 'mixed'")

    self.setParameter("printfe","mbar")
    self.setParameter("dualtopology1","1 2 synctrans syncrot")

    if softcore != 'none':
      softcore_options1, softcore_options2 = "solute 1", "solute 2"
      if softcore == 'mixed':
        softcore_options1 += " atoms "
        softcore_options2 += " atoms "

        tem_file = sim.TemplateFile(filename=templates[0])
        cmap = {}
        make_single._auto_map(tem_file.templates[0], tem_file.templates[1],
                              sim.PDBFile(filename=solutes[0]),
                              sim.PDBFile(filename=solutes[1]),
                              cmap)

        for i, atom in enumerate(tem_file.templates[0].atoms):
          if atom.name not in cmap:
            softcore_options1 += "%d " % (i+1)

        for i, atom in enumerate(tem_file.templates[1].atoms):
          if atom.name not in cmap.values():
            softcore_options2 += "%d " % (i+1)

      self.setParameter("softcore1", softcore_options1)
      self.setParameter("softcore2", softcore_options2)
      self.setParameter("softcoreparams","coul 1 delta 0.2 deltacoul 2.0 power 6 soft66")
    self.setParameter("dlambda","0.001")
    self.setParameter("lambdare","%d %s"%(2*dumpfreq," ".join("%.3f"%l for l in lambdaval)))

    self.setDump("results write results",dumpfreq)
    self.setDump("results writeinst results_inst",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)

    if protein is not None :
      for restsol in restrained :
        pdbobj = sim.PDBFile(filename=solutes[restsol])
        for ind,tem in enumerate(templates) :
          temobj = sim.TemplateFile(filename=templates[ind])
          for mol_template in temobj.templates :
            if mol_template.name in pdbobj.header :
              resname = pdbobj.residues[1].name
              resatom = str(mol_template.atoms[0]).strip().split()
              for atom in pdbobj.residues[1].atoms :
                if atom.name in resatom[1] : atmcoords = atom.coords
        self.setChunk("id add %d solute %d %s %s"%(restsol+1,restsol+1,resatom[1],resname))
        self.setChunk("restraint add %d cartesian harmonic %.3f %.3f %.3f %d"%(restsol+1,atmcoords[0],atmcoords[1],atmcoords[2],FORCE_CONSTANT))
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"standard",self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class SingleTopology(ProteinLigandSimulation) :
  """ 
  This a command file for a single-topology protein-ligand simulation
  Generates input for proteins, solutes and solvent
  write lambda-replica exchange,
  dump results, pdb and restart files
  equilibrate and production run
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute1.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    nequil=5E6,
                    nprod=40E6,
                    dumpfreq=1E5,
                    lambdaval=None,
                    ranseed=None,
                    outfolder="out",
                    restrained=[]) :  
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    lambdaval : float or list of floats
      the lambda values to perform the simulation at
    outfolder  : string, optional
      the folder for all output files
    restrained : string, optional
      the solutes on which restrains should be applied
    
    Raises
    ------
    SetupError
      no lambda values given
    """      
    if len(outfolder) == 0 : outfolder = "out"          
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)
     
    if lambdaval is None or len(lambdaval) < 2 :
      raise sim.SetupError("Must give at least two lambda values")

    self.setParameter("printfe","mbar")
    self.setParameter("dlambda","0.001")
    self.setParameter("lambdare","%d %s"%(2*dumpfreq," ".join("%.3f"%l for l in lambdaval)))

    self.setDump("results write results",dumpfreq)
    self.setDump("results writeinst results_inst",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"standard",self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class RestraintRelease(ProteinLigandSimulation) :
  """ 
  Command file to calculate the free energy of introducing a harmonic restraint
  during the bound leg of a simulation.
  """
  
  def __init__(self,protein="protein.pdb",
                    solutes=["solute1.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    nequil=5E6,
                    nprod=40E6,
                    dumpfreq=1E5,
                    lambdaval=None,
                    ranseed=None,
                    outfolder="out",
                    restrained=[]) :
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    lambdaval : float or list of floats
      the lambda values to perform the simulation at
    outfolder  : string, optional
      the folder for all output files
    restrained : string, optional
      the solutes on which restrains should be applied
    
    Raises
    ------
    SetupError
      no lambda values given
    """
    if len(outfolder) == 0 : outfolder = "out"          
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)
     
    if lambdaval is None or len(lambdaval) < 2 :
      raise sim.SetupError("Must give at least two lambda values")

    self.setParameter("printfe","mbar")
    self.setParameter("dlambda","0.001")
    self.setParameter("lambdare","%d %s"%(2*dumpfreq," ".join("%.3f"%l for l in lambdaval)))

    self.setDump("results write results",dumpfreq)
    self.setDump("results writeinst results_inst",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    for restsol in restrained :
      pdbobj = sim.PDBFile(filename=solutes[restsol])
      for ind,tem in enumerate(templates) :
        temobj = sim.TemplateFile(filename=templates[ind])
        for mol_template in temobj.templates :
          if mol_template.name in pdbobj.header :
            resname = pdbobj.residues[1].name
            resatom = str(mol_template.atoms[0]).strip().split()
            for atom in pdbobj.residues[1].atoms :
              if atom.name in resatom[1] : atmcoords = atom.coords
      self.setChunk("id add %d solute %d %s %s"%(restsol+1,restsol+1,resatom[1],resname))
      self.setChunk("restraint add %d cartesian harmonic %.3f %.3f %.3f %d lambda"%(restsol+1,atmcoords[0],atmcoords[1],atmcoords[2],FORCE_CONSTANT))
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"standard",self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

    
#
# Helper routine
#

def _setbox(simulation,waters,inbox) :

  waters = sim.PDBFile(filename=waters)
  headerbox = sim.find_box(waters)
      
  if headerbox is None :
    if inbox is None :
      raise sim.SetupError("Cannot setup simulation without a box")
    pdbobj = sim.PDBFile(filename=inbox) 
    outbox = pdbobj.getBox()    
  else :
    outbox = headerbox

  if "origin" in outbox :
    for i,param in enumerate(["x","y","z"]) :
      simulation.setParameter("origin"+param,outbox["origin"][i])
  elif "center" in outbox :
    for i,param in enumerate(["x","y","z"]) :
      simulation.setParameter("center"+param,outbox["center"][i])
  for i,param in enumerate(["x","y","z"]) :
    simulation.setParameter(param,outbox["len"][i])

class GCMC(ProteinLigandSimulation) :
  """ 
  This a command file for a GCMC simulation
  Generates input for proteins, solutes and solvent
  write GCMC parameters
  dump results, pdb and restart files
  equilibrate and production run
  """
  def __init__(self,protein="protein.pdb",
                    solutes=[],
                    solvent="water.pdb",
                    templates=[],
                    gcmcwater="gcmc_water.pdb",
                    gcmcbox=None,   
                    nequil=5E6,
                    nprod=80E6,
                    dumpfreq=2E5,
                    adamval=None,
                    ranseed=None,
                    watmodel="tip4p",
                    outfolder="out_gcmc") :  
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    gcmcwater : string, optional
      filename of the water to do gcmc one
    gcmcbox : dictionary of numpy array, optional
      defines the box if not found in gcmcwater
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    adamval : list of floats
      the Adams values to perform the simulation at
    outfolder  : string, optional
      the folder for all output files
    
    Raises
    ------
    SetupError
      no protein given
      no box dimensions found
      no GCMC water given
      no Adam values givens
    """      
    #if len(outfolder) == 0 : outfolder = "out_gcmc"
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)
     
    if adamval is None :
      raise sim.SetupError("Must give at least one Adam value")

    if protein is None :
      raise sim.SetupError("Cannot setup GCMC without protein")

    if gcmcwater is None  :
      raise sim.SetupError("Cannot setup GCMC without any GCMC water")
 
    self.setParameter("#"," GCMC specific parameters")
    self.setParameter("gcmc","0")
    if watmodel == "tip4p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip4p.tem"))
    elif watmodel == "tip3p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip3p.tem"))
    self.setParameter("grand1",gcmcwater)
    #self.setParameter("outfolder",outfolder)
    if isinstance(adamval,float) or isinstance(adamval,int):
      self.setParameter("potential","%.3f"%adamval)
    elif len(adamval) == 1 :
      self.setParameter("potential","%.3f"%adamval[0])
      #self.setParameter("outfolder",outfolder)
    else :
      self.setParameter("multigcmc","%d %s " % (2*dumpfreq," ".join("%.3f"%a for a in adamval) ))
      #self.setParameter("refolder",outfolder)
    _setbox(self,gcmcwater,gcmcbox)
    self.setParameter("#"," End of GCMC specific parameters")
  
    self.setDump("results write results",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"gcmc",self.periodic)
    self.setChunk("equilibrate %d solvent=0 protein=0 solute=0 insertion=333 deletion=333 gcsolute=333" % nequil)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class Jaws1(ProteinLigandSimulation) :
  """ 
  This a command file for a JAWS-1 simulation
  Generates input for proteins, solutes and solvent
  write JAWS-1 parameters
  dump results, pdb and restart files
  equilibrate and production run
  """
  def __init__(self,protein="protein.pdb",
                    solutes=[],
                    solvent="water.pdb",
                    templates=[],
                    outfolder="",
                    ranseed=None,
                    jawswater="jaws_water.pdb",
                    jawsbox=None,   
                    nequil=5E6,
                    nprod=40E6,
                    watmodel="tip4p",
                    dumpfreq=1E5) :  
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    jawswater : string, optional
      filename of the water to do JAWS-1 one
    jawsbox : dictionary of numpy array, optional
      defines the box if not found in jawswater
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    outfolder  : string, optional
      the folder for all output files
    
    Raises
    ------
    SetupError
      no protein given
      no box dimensions found
      no JAWS water given
    """      
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)
     
    if protein is None :
      raise sim.SetupError("Cannot setup GCMC without protein")

    if jawswater is None  :
      raise sim.SetupError("Cannot setup JAWS-1 without any JAWS water")
 
    self.setParameter("#"," JAWS-1 specific parameters")
    self.setParameter("jaws1","0")
    if watmodel == "tip4p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip4p.tem"))
    elif watmodel == "tip3p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip3p.tem"))
    self.setParameter("grand1",jawswater)
    _setbox(self,jawswater,jawsbox)
    self.setParameter("#"," End of JAWS specific parameters")
  
    self.setDump("results write results",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"jaws1",self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class Jaws2(ProteinLigandSimulation) :
  """ 
  This a command file for a JAWS-2 simulation
  Generates input for proteins, solutes and solvent
  write JAWS-2 parameters
  dump results, pdb and restart files
  equilibrate and production run
  """
  def __init__(self,protein="protein.pdb",
                    solutes=[],
                    solvent="water.pdb",
                    templates=[],
                    outfolder="",
                    ranseed=None,
                    jawswater="jaws_water.pdb",
                    jbias=6.5,
                    nequil=5E6,
                    nprod=40E6,
                    watmodel="tip4p",
                    dumpfreq=1E5) :  
    """
    Parameters
    ----------
    protein : string, optional
      the filename of a protein pdb file
    solutes : list of strings, optional
      filenames of solute pdb files
    solvent : string, optional
      the filename of a solvent pdb file
    templates : list of strings, optional
      filenames of template files to be included
    outfolder  : string, optional
      the folder for all output files
    jawswater : string, optional
      filename of the water to do JAWS-2 one
    jbias : float, optional
      the bias
    nequil : int, optional
      number of equilibration moves
    nprod : int, optional
      number of production moves
    dumfreq : int, optional
      the dump frequency 
    
    Raises
    ------
    SetupError
      no protein given
      no box dimensions found
      no JAWS water given
    """      
    solvent,jawssolvent = solvent.split()
    if len(outfolder) == 0 : outfolder = "out_jaws2"
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates,outfolder=outfolder,ranseed=ranseed)
     
    if protein is None :
      raise sim.SetupError("Cannot setup GCMC without protein")

    if jawswater is None  :
      raise sim.SetupError("Cannot setup JAWS-2 without any JAWS water")
 
    self.setParameter("#"," JAWS-2 specific parameters")
    self.setParameter("jaws2","1")
    if isinstance(jbias,float) or isinstance(jbias,int):
      self.setParameter("jbias","%.3f"%jbias)
    elif len(jbias) == 1 :
      self.setParameter("jbias","%.3f"%jbias[0])
    else :
      self.setParameter("multijaws2"," ".join("%.3f"%j for j in jbias))
    if watmodel == "tip4p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip4p.tem"))
    elif watmodel == "tip3p" :
      self.setForceField(os.path.join("$PROTOMSHOME","data","gcmc_tip3p.tem"))
    self.setSolvent(2,jawssolvent)
    self.setParameter("grand1",jawswater)
    watobj = sim.PDBFile(filename=jawswater)
    ocoord = watobj.residues[1].atoms[0].coords
    for i,param in enumerate(["x","y","z"]) :
      self.setParameter("origin"+param,ocoord[i]-1.5)
    for i,param in enumerate(["x","y","z"]) :
      self.setParameter(param,"3.0")
    self.setParameter("#"," End of JAWS specific parameters")
  
    self.setDump("results write results",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = _assignMoveProbabilities(protein,solutes,solvent,"jaws2",self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

def generate_input(protein,ligands,templates,protein_water,ligand_water,ranseed,settings) :

  """
  Generates ProtoMS command files
  
  If protein is a valid filename a protein(-ligand) simulation is setup, and if
  ligands contains at least one valid filename a ligand simulation is setup. 
  
  Parameters
  ----------
  protein : string
    the filename of a protein pdb file
  ligands : list of strings
    filenames of solute pdb files
  templates : list of strings
    filenames of template files to be included
  protein_water : string
    the filename of a solvent pdb file for the protein
  ligand_water : string
    the filename of a solvent pdb file for the ligands
  settings : Namespace object (from argparse)
    additional settings

  Returns
  -------
  free_cmd : ProtoMSSimulation
    the command file for the ligand simulation
  bnd_cmd : ProtoMSSimulation
    the command file for the protein simulation
  gas_cmd : ProtoMSSimulation
    the command file for the gas simulation
  """
  
  logger.debug("Running generate_input with arguments: ")
  logger.debug("\tprotein       = %s"%protein) 
  if ligands is not None :
    logger.debug("\tligands       = %s"%" ".join(ligands)) 
  else :
    logger.debug("\tligand        = None")
  if templates is not None :
    logger.debug("\ttemplates     = %s"%" ".join(templates))
  else :
    logger.debug("\ttemplates     = None")
  logger.debug("\tprotein_water = %s"%protein_water)
  logger.debug("\tligand_water  = %s"%ligand_water)
  logger.debug("\tsettings      = %s"%settings)
  logger.debug("This will make an input file for ProtoMS")
  
  free_cmd = bnd_cmd = gas_cmd = None
  ranseed=settings.ranseed
  ProteinLigandSimulation.gaffversion = settings.gaff

  if settings.simulation == "equilibration" :

    if protein is not None : # Setup a protein simulation
      bnd_cmd = Equilibration(protein=protein,solutes=ligands, 
                              templates=templates,solvent=protein_water,outfolder=settings.outfolder+"_bnd",
                              nsteps=settings.nequil,pdbfile="equil_bnd.pdb",ranseed=ranseed)
    elif ligands is not None :      
      if settings.dovacuum : 
        solvent = None
        prestr = "gas"
      else :
        solvent = ligand_water
        prestr = "free"

      free_cmd = Equilibration(protein=None,solutes=ligands[:min(len(ligands),2)], 
                              templates=templates,
                              solvent=solvent,outfolder=settings.outfolder+"_"+prestr,
                              nsteps=settings.nequil,pdbfile="equil_%s.pdb"%prestr,ranseed=ranseed)
      if settings.dovacuum :
        gas_cmd = free_cmd
        free_cmd = None
      
  elif settings.simulation == "sampling" :
    if protein is not None : # Setup a protein simulation
      bnd_cmd = Sampling(protein=protein,solutes=ligands, 
                         templates=templates,solvent=protein_water,
                         outprefix="bnd_",
                         nequil=settings.nequil,outfolder=settings.outfolder+"_bnd",
                         nprod=settings.nprod,dumpfreq=settings.dumpfreq,ranseed=ranseed,tune=settings.tune)
    elif ligands is not None :      
      if settings.dovacuum : 
        solvent = None
        prestr = "gas"
      else :
        solvent = ligand_water
        prestr = "free"

      free_cmd = Sampling(protein=None,solutes=ligands[:min(len(ligands),2)], 
                         templates=templates,
                         solvent=solvent,outprefix=prestr+"_",
                         nequil=settings.nequil,outfolder=settings.outfolder+"_"+prestr,
                         nprod=settings.nprod,dumpfreq=settings.dumpfreq,ranseed=ranseed,tune=settings.tune)
      if settings.dovacuum :
        gas_cmd = free_cmd
        free_cmd = None
          
  elif settings.simulation in ["dualtopology","singletopology"] :
  
    cmdcls = {"dualtopology":DualTopology,"singletopology":SingleTopology}
    extra_kwargs = {}
    if settings.simulation == "dualtopology":
      extra_kwargs['softcore'] = settings.softcore
      
    if ligands is None :
      raise tools.SetupError("No ligands loaded, cannot do dual-topology simulations.")

    if len(settings.lambdas) == 1 :
      nlambdas = int(settings.lambdas[0])
      lambdavals = np.linspace(0,1,nlambdas,endpoint=True)
    else :
      lambdavals = settings.lambdas
      nlambdas = len(lambdavals)
    logger.info("")
    logger.info("Will simulate with %s lambda values"%nlambdas)
  
    if hasattr(settings,"outfolder") and settings.outfolder != "" :
      outfolder = settings.outfolder
    else :
      outfolder = "out"

    rest_solutes = []
    if settings.simulation == 'dualtopology' and settings.absolute:
      rest_solutes.append(0) 
      rest_solutes.append(1)

    if settings.simulation == "singletopology"  :
      gas_cmd = cmdcls[settings.simulation](protein=None,solutes=ligands[:min(len(ligands),2)], 
                              templates=templates,solvent=None,
                              lambdaval=lambdavals,nequil=settings.nequil,ranseed=ranseed,
                                            nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder+"_gas",restrained=rest_solutes_kwargs)
    elif settings.simulation == "dualtopology" and settings.absolute:
      #make use of the unused gas_cmd to produce a cmd file to calculate the effect of introducing a harmonic restraint
      gas_cmd = RestraintRelease(protein=protein,solutes=ligands[:1],templates=templates[:1],solvent=protein_water,
                                 lambdaval=(0.,1.),nequil=settings.nequil,ranseed=ranseed,
                                 nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder+"_bnd_rstr",restrained=rest_solutes[:1])

    free_cmd = cmdcls[settings.simulation](protein=None,solutes=ligands[:min(len(ligands),2)], 
                            templates=templates,solvent=ligand_water,
                            lambdaval=lambdavals,nequil=settings.nequil,ranseed=ranseed,
                            nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder+"_free",restrained=rest_solutes,**extra_kwargs)
      
    if protein is not None :
      bnd_cmd = cmdcls[settings.simulation](protein=protein,solutes=ligands, 
                             templates=templates,solvent=protein_water,
                             lambdaval=lambdavals,nequil=settings.nequil,ranseed=ranseed,
                             nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder+"_bnd",restrained=rest_solutes,**extra_kwargs)

  elif settings.simulation == "gcmc" :
  
     
      if hasattr(settings,"outfolder") and settings.outfolder != "" :
        outfolder = settings.outfolder
      else :
        outfolder = "out_gcmc"

      if settings.adamsrange is not None:
        if len(settings.adamsrange) not in ( 2, 3 ):
          msg = "If using --adamsrange, please specify exactly two Adams values and optionally the number of desired values "
          logger.error(msg)
          raise sim.SetupError(msg)
        else:
          adams1 = settings.adamsrange[0]
          adams2 = settings.adamsrange[1]
          try:
            nadams = settings.adamsrange[2]
          except IndexError:
            nadams = int ( abs ( adams1 - adams2 ) ) + 1
          if adams2-adams1 > 0 : 
            settings.adams = np.linspace(adams2,adams1,nadams)
          else: 
            settings.adams = np.linspace(adams1,adams2,nadams)

      bnd_cmd = GCMC(protein=protein,solutes=ligands, 
                     templates=templates,solvent=protein_water,gcmcwater=settings.gcmcwater,
                     adamval=settings.adams,nequil=settings.nequil,gcmcbox=settings.gcmcbox,
                     nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder,ranseed=ranseed,watmodel=settings.watmodel)    

  elif settings.simulation == "jaws1" :
  
      bnd_cmd = Jaws1(protein=protein,solutes=ligands, 
                     templates=templates,solvent=protein_water,jawswater=settings.gcmcwater,
                     nequil=settings.nequil,jawsbox=settings.gcmcbox,ranseed=ranseed,
                     nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=settings.outfolder) 

  elif settings.simulation == "jaws2" :
  
      if hasattr(settings,"outfolder") and settings.outfolder != "" :
        outfolder = settings.outfolder
      else :
        outfolder = "out_jaws2"

      bnd_cmd = Jaws2(protein=protein,solutes=ligands, 
                     templates=templates,solvent=protein_water,jawswater=settings.gcmcwater,
                     nequil=settings.nequil,jbias=settings.jawsbias,ranseed=ranseed,
                     nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder=outfolder) 
                     
  return free_cmd,bnd_cmd,gas_cmd  

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to create a ProtoMS command file")
  parser.add_argument('-s','--simulation',choices=["sampling","equilibration","dualtopology","singletopology","gcmc","jaws1","jaws2"],help="the kind of simulation to setup",default="equilibration")
  parser.add_argument('--dovacuum',action='store_true',help="turn on vacuum simulation for simulation types equilibration and sampling",default=False)
  parser.add_argument('-p','--protein',help="the name of the protein file")
  parser.add_argument('-l','--ligands',nargs="+",help="the name of the ligand pdb files")
  parser.add_argument('-t','--templates',nargs="+",help="the name of ProtoMS template files")
  parser.add_argument('-pw','--protwater',help="the name of the solvent for protein")
  parser.add_argument('-lw','--ligwater',help="the name of the solvent for ligand")
  parser.add_argument('-o','--out',help="the prefix of the name of the command file",default="run")
  parser.add_argument('--outfolder',help="the ProtoMS output folder",default="out")
  parser.add_argument('--gaff',help="the version of GAFF to use for ligand",default="gaff16")
  parser.add_argument('--lambdas',nargs="+",type=float,help="the lambda values or the number of lambdas",default=[16])
  parser.add_argument('--adams',nargs="+",type=float,help="the Adam/B values for the GCMC",default=0)
  parser.add_argument('--adamsrange',nargs="+",type=float,help="the upper and lower Adam/B values for the GCMC and, optionally, the number of values desired (default value every 1.0), e.g. -1 -16 gives all integers between and including -1 and -16",default=None)
  parser.add_argument('--jawsbias',nargs="+",type=float,help="the bias for the JAWS-2",default=0)
  parser.add_argument('--gcmcwater',help="a pdb file with a box of water to do GCMC on")
  parser.add_argument('--gcmcbox',help="a pdb file with box dimensions for the GCMC box")
  parser.add_argument('--watmodel',help="the name of the water model. Default = tip4p",choices=[ 'tip3p', 'tip4p'],default='tip4p')
  parser.add_argument('--nequil',type=float,help="the number of equilibration steps",default=5E6)
  parser.add_argument('--nprod',type=float,help="the number of production steps",default=40E6)
  parser.add_argument('--dumpfreq',type=float,help="the output dump frequency",default=1E5)
  parser.add_argument('--absolute',action='store_true',help="whether an absolute free energy calculation is to be run. Default=False",default=False)
  parser.add_argument('--ranseed',help="the value of the random seed you wish to simulate with. If None, then a seed is randomly generated. Default=None",default=None)
  parser.add_argument('--tune',action='store_true',help=argparse.SUPPRESS,default=False)
  parser.add_argument('--softcore', type=str, default='mixed', 
                      choices=('mixed', 'all', 'none'),
                      help="determine which atoms to apply softcore potentials to.\n "
                           "'all'=softcores applied to all atoms of both solutes, "
                           "'none'=softcores not applied to any atoms\n "
                           "'mixed'=softcores will be applied only to non matching "
                           "atoms within ligand structures")
  return parser

if __name__ == "__main__":

  args = get_arg_parser().parse_args()

  # Setup the logger
  logger = sim.setup_logger("generate_input_py.log")
  
  free_cmd,bnd_cmd,gas_cmd = generate_input(args.protein,args.ligands,args.templates,args.protwater,args.ligwater,args.ranseed,args) 

  args.out = args.out.lower() #protoMS cannot handle cmd files containing upper case letters
  if free_cmd is not None : 
    free_cmd.writeCommandFile(args.out+"_free.cmd")
  if bnd_cmd is not None : 
    if args.simulation == "gcmc" :
      bnd_cmd.writeCommandFile(args.out+"_gcmc.cmd") 
    elif args.simulation in ["jaws1","jaws2"] :
      bnd_cmd.writeCommandFile(args.out+"_jaws.cmd") 
    else :
      bnd_cmd.writeCommandFile(args.out+"_bnd.cmd")      
  if gas_cmd is not None : 
    gas_cmd.writeCommandFile(args.out+"_gas.cmd")    
