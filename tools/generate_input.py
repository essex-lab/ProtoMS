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

import numpy as np

import simulationobjects


def assignMoveProbabilities(protein,solute,solvent,isgcmc,isperiodic) :
  """ Assigns move probabilities for protein, solute, solvent and box
      Does this by some "heuristic" rule
  """
  pvolume = 0.0
  numprot = 0.0
  numsolu = 0.0
  numsolv = 0.0
  pinsert = 0.0
  pdelete = 0.0
  pgcsolu = 0.0

  accel_prot = 5					# The amount by which the protein be accelerated over the solvent.
  accel_solu = accel_prot*10		# The amount by which the solute be accelerated over the solvent.
  addto = 1000.0					# Move proportion add up to this number. This is just for clarity.

  if protein is not None:
    protobj = simulationobjects.PDBFile(filename=protein)
    numprot = float(len(protobj.residues))*accel_prot
  if solute is not None:
    soluobj =  simulationobjects.PDBFile(filename=solute[0])
    numsolu = float(len(soluobj.residues))*accel_solu
  if solvent is not None:
    solvobj =  simulationobjects.PDBFile(filename=solvent)
    numsolv = float(len(solvobj.solvents))
  
  numtot = numprot + numsolu + numsolv
  if isperiodic:
    pvolume = int(0.002*addto) + 1
  else:
    pvolume = 0
  if isgcmc == False:
    psolvent = int(round(addto*numsolv/numtot)) - pvolume
    psolute = int(round(addto*numsolu/numtot))
    pprotein = int(round(addto*numprot/numtot))
    return "solvent=%d protein=%d solute=%d volume=%d"%(psolvent,pprotein,psolute,pvolume)
  else:
    pinsert = int(round(addto/6))				# Giving half of all moves to GCMC moves.
    pdelete = int(round(addto/6))
    pgcsolu = int(round(addto/6))
    psolvent = int(round(addto*numsolv/numtot/2))
    psolute = int(round(addto*numsolu/numtot/2))
    pprotein = int(round(addto*numprot/numtot/2))
    return "solvent=%d protein=%d solute=%d insertion=%d deletion=%d gcsolute=%d"%(psolvent,pprotein,psolute,pinsert,pdelete,pgcsolu)

def assignMoveProbabilities_old(hasprotein,hassolute,hassolvent,isperiodic) :
  """ Assigns move probabilities for protein, solute, solvent and box
      Does this by some "heuristic" rule
  """
  psolute = 0
  psolvent = 0
  pprotein = 0
  pvolume = 0
  if hassolvent :
    if hasprotein :
      psolvent= 857
    else  :
      psolvent = 990
    if isperiodic : 
      pvolume = 1
  if hasprotein :
    pprotein = 128
  if hassolute :
    if hasprotein :
      psolute = 14
    else :
      psolute = 9
  return "solvent=%d protein=%d solute=%d volume=%d"%(psolvent,pprotein,psolute,pvolume)

class ProtoMSSimulation :
    """This is a ProtoMS command file object.This object holds all of the information
       about the simulation and can either run the simulation with this
       information, or can write this information to a command file.
       
    Attributes:
        lines -- all the lines of the ProtoMS command file    
    """
    def __init__(self,filename=None):
        self._lines = []
        if (not filename is None):
            self.readCommandFile(filename)        
    def __str__(self):
        string = "\n".join(self._lines)           
        return string        
    def _setkey(self,key,value):
      self._lines.append("%s %s"%(key,value))        
    def clear(self):
        """Clear the simulation description"""    
        self._lines = []
    def readCommandFile(self,filename):
        """Read in a simulation from a command file"""
        lines = open(filename,"r").readlines()
     
        # Attempts to read PROTOMSHOME environmental variable
        pmhome = os.getenv("PROTOMSHOME")       
 
        for line in lines:
            #split the line into words
            words = line.split()
            #skip short lines or comments
            if (len(words) < 2): continue
            if (words[0].find("#") == 0): continue

            arguments = ' '.join(words[1:])
            if pmhome != None :
              arguments = arguments.replace("$PROTOMSHOME",pmhome)
            self._setkey(words[0],arguments)
    def setStream(self,stream,filename):
        """Set the direction of the output stream to a filename"""
        key = "stream%s" % stream.lower()
        self._setkey(key,filename)
    def setParameter(self,parameter,value):
        """Set the parameter to value"""
        self._setkey(parameter.lower(),value)
    def setChunk(self,chunk):
        """This adds the simulation chunk to 'chunk'. """
        self._setkey("chunk",chunk)
    def setDump(self,dump,freq):
        """This adds the simulation dump to 'dump' with frequnecy freq"""
        self._setkey("dump","%d %s"%(freq,dump))
    def setProtein(self,n,filename):
        """Set the filename for proteinN. It then returns n+1"""
        key = "protein%d" % n
        self._setkey(key,filename)
        return n+1      
    def setSolute(self,n,filename):
        """Set the filename for soluteN. It then returns n+1"""
        key = "solute%d" % n
        self._setkey(key,filename)
        return n+1
    def setSolvent(self,n,filename):
        """Set the filename for solventN. It then returns n+1"""
        key = "solvent%d" % n
        self._setkey(key,filename)
        return n+1
    def setForceField(self,filename):
        """Set the filename for parfile"""
        self._setkey("parfile",filename)
    def writeCommandFile(self,filename):
        """Write the contents of this simulation object to the
           command file 'filename'"""
        f = open(filename,"w")
        for line in self._lines:
            f.write("%s\n" % line)            
        f.close()
    def run(self,exe,cmdfile=None):
        """Run this simulation using the executable 'exe' - returns the exit status
           of the run"""
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
  """ This a generic command file for a protein-ligand simulation
      Generates input for proteins, solutes and solvent
      but does not write any chunks
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"]) :
    ProtoMSSimulation.__init__(self)
    self.setForceField("$PROTOMSHOME/parameter/amber99.ff")
    self.setForceField("$PROTOMSHOME/parameter/solvents.ff")
    self.setForceField("$PROTOMSHOME/parameter/amber99-residues.ff")
    self.setForceField("$PROTOMSHOME/parameter/gaff.ff")
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
    self.setStream("header","off")
    self.setStream("detail","off")
    for s in ["warning","info","fatal","results"] :
      self.setStream(s,s)
    self.setParameter("cutoff",10.0)
    self.setParameter("feather",0.5)
    if not solvent is None :
      self.setParameter("boundary","solvent")
      if self.periodic  :
        self.setParameter("pressure","1")

class Equilibration(ProteinLigandSimulation) :
  """ This a command file for equilibration of a protein-ligand system
      Generates input for proteins, solutes and solvent
      and write an equilibration and a pdb chunk
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    nsteps=1E5,
                    pdbfile="equilibrated.pdb") :
                    
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates)

    #moves = assignMoveProbabilities(protein is not None,solutes,solvent is not None,self.periodic)
    moves = assignMoveProbabilities(protein,solutes,solvent,False,self.periodic)
    self.setChunk("equilibrate %d %s"%(nsteps,moves))
    self.setChunk("pdb all solvent=all file=%s standard"%pdbfile)

class Sampling(ProteinLigandSimulation) :
  """ This a command file for MC sampling of a protein-ligand system
      Generates input for proteins, solutes and solvent
      and write an equilibration and a simulate
  """
  def __init__(self,protein="protein.pdb",
                    solutes=["solute.pdb"],
                    solvent="water.pdb",
                    templates=["solute.tem"],
                    nequil=5E6,
                    nprod=40E6,
                    dumpfreq=1E5,
                    outprefix="") :
                    
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates)

    self.setDump("results write %sresults"%outprefix,dumpfreq)
    self.setDump("pdb all solvent=all file=%sall.pdb standard"%outprefix,dumpfreq)
    self.setDump("restart write %srestart"%outprefix,dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = assignMoveProbabilities(protein,solutes,solvent,False,self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

class DualTopology(ProteinLigandSimulation) :
  """ This a command file for a dual-topology protein-ligand simulation
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
                    lambdaval=None,
                    outfolder="out") :  
                    
    ProteinLigandSimulation.__init__(self,protein=protein,solutes=solutes,solvent=solvent,templates=templates)

    if len(solutes) < 2 :
      raise simulationobjects.SetupError("Cannot do dual topology with less than 2 solutes")
      
    if lambdaval is None or len(lambdaval) < 2 :
      raise simulationobjects.SetupError("Must give at least two lambda values")

    self.setParameter("pdbparams","on")
    self.setParameter("printfe","mbar")
    self.setParameter("dualtopology1","1 2 synctrans syncrot")
    self.setParameter("softcore1","solute 1")
    self.setParameter("softcore2","solute 2")
    self.setParameter("softcoreparams","coul 1 delta 0.2 deltacoul 2.0 power 6")
    self.setParameter("lambda","0.500 0.501 0.499")
    self.setParameter("lambdare","%d %s"%(2*dumpfreq," ".join("%.3f"%l for l in lambdaval)))
    self.setParameter("refolder",outfolder)

    self.setDump("results write results",dumpfreq)
    self.setDump("pdb all solvent=all file=all.pdb standard",dumpfreq)
    self.setDump("restart write restart",dumpfreq)
    self.setDump("averages reset",dumpfreq)
    
    moves = assignMoveProbabilities(protein,solutes,solvent,False,self.periodic)
    self.setChunk("equilibrate %d %s"%(nequil,moves))        
    self.setChunk("simulate %d %s"%(nprod,moves))

def generate_input(protein,ligands,templates,protein_water,ligand_water,settings) :

  free_cmd = bnd_cmd = False
  
  if settings.simulation == "equilibration" :
    if ligands is not None :
      if settings.dovacuum : 
        solvent = None
      else :
        solvent = ligand_water
      free_cmd = Equilibration(protein=None,solutes=ligands[:min(len(ligands),2)], 
                              templates=templates,
                              solvent=solvent,
                              nsteps=settings.nequil,pdbfile="equil_free.pdb")
      
    if protein is not None :
      if settings.dovacuum :
        solvent = None
      else :
        solvent = protein_water  
      bnd_cmd = Equilibration(protein=protein,solutes=ligands, 
                              templates=templates,solvent=solvent,
                              nsteps=args.nequil,pdbfile="equil_bnd.pdb")
      
  elif settings.simulation == "sampling" :
    if ligands is not None :
      if args.dovacuum : 
        solvent = None
      else :
        solvent = ligand_water
      free_cmd = Sampling(protein=None,solutes=ligands[:min(len(ligands),2)], 
                         templates=templates,
                         solvent=solvent,outprefix="free_",
                         nequil=settings.nequil,
                         nprod=settings.nprod,dumpfreq=settings.dumpfreq)
      
    if protein is not None :
      if settings.dovacuum :
        solvent = None
      else :
        solvent = protein_water  
      bnd_cmd = Sampling(protein=protein,solutes=ligands, 
                         templates=templates,solvent=solvent,
                         outprefix="bnd_",
                         nequil=settings.nequil,
                         nprod=settings.nprod,dumpfreq=settings.dumpfreq)
            
  elif settings.simulation == "dualtopology" :
  
    if ligands is None :
      raise tools.SetupError("No ligands loaded, cannot do dual-topology simulations.")

    if len(settings.lambdas) == 1 :
      nlambdas = int(settings.lambdas[0])
      lambdavals = np.linspace(0,1,nlambdas,endpoint=True)
    else :
      lambdavals = settings.lambdas
      nlambdas = len(lambdavals)
    print "\nWill simulate with %s lambda values"%nlambdas
  
    free_cmd = DualTopology(protein=None,solutes=ligands[:min(len(ligands),2)], 
                            templates=templates,solvent=ligand_water,
                            lambdaval=lambdavals,nequil=settings.nequil,
                            nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder="out_free")
      
    if protein is not None :
      bnd_cmd = DualTopology(protein=protein,solutes=ligands, 
                             templates=templates,solvent=protein_water,
                             lambdaval=lambdavals,nequil=settings.nequil,
                             nprod=settings.nprod,dumpfreq=settings.dumpfreq,outfolder="out_bnd")
                        
  return free_cmd,bnd_cmd  

if __name__ == "__main__":

  import argparse
  import numpy as np

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to create a ProtoMS command file")
  parser.add_argument('-s','--simulation',choices=["vacuum","equilibration","dualtopology"],help="the kind of simulation to setup",default="equilibration")
  parser.add_argument('-p','--protein',help="the name of the protein file")
  parser.add_argument('-l','--ligands',nargs="+",help="the name of the ligand pdb files")
  parser.add_argument('-t','--templates',nargs="+",help="the name of ProtoMS template files")
  parser.add_argument('-w','--water',help="the name of the solvent")
  parser.add_argument('-o','--out',help="the name of the command file",default="run.cmd")
  parser.add_argument('--lambdas',nargs="+",type=float,help="the lambda values or the number of lambdas",default=[16])
  parser.add_argument('--nequil',type=int,help="the number of equilibration steps",default=5E6)
  parser.add_argument('--nprod',type=int,help="the number of production steps",default=40E6)
  parser.add_argument('--dumpfreq',type=int,help="the output dump frequency",default=1E5)
  args = parser.parse_args()

  if args.simulation == "vacuum" :
    cmdfile = Vacuum(protein=args.protein,solutes=args.ligands,templates=args.templates,nsteps=args.nprod)  
  elif args.simulation == "equilibration" :
    cmdfile = Equilibration(protein=args.protein,solutes=args.ligands,solvent=args.water,templates=args.templates,nsteps=args.nequil)
  elif args.simulation == "dualtopology" : 
    if len(args.ligands) < 2 :
      simulationobjects.SetupError("Cannot perform dual topology with less than two ligands.")
    if len(args.lambdas) == 1 :
      nlambdas = int(args.lambdas[0])
      lambdavals = np.linspace(0,1,nlambdas,endpoint=True)
    else :
      lambdavals = args.lambdas
      nlambdas = len(lambdavals)
    cmdfile = DualTopology(protein=args.protein,solutes=args.ligands,solvent=args.water,templates=args.templates,nequil=args.nequil,nprod=args.nprod,dumpfreq=args.dumpfreq,lambdaval=lambdavals)
  cmdfile.writeCommandFile(args.out)
