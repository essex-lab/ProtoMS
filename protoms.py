# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Program to setup, run and analyse a ProtoMS simulation
"""

import argparse
import os
import subprocess

import numpy as np

import tools
from tools import simulationobjects

def get_prefix(filename) :
  h,t = os.path.splitext(filename)
  return h
  
def locate_file(filename,folders) :
  """ Find a file 
  """
  # Try to see if the filename as given exists
  if os.path.isfile(filename) :
    return filename
  else :
    # Remove everything but the actual filename
    h,t = os.path.split(filename)
    # Loop over all folder and try to find it there
    for f in folders :
      test = os.path.join(f,t)
      if os.path.isfile(test) :
        return test
  # If we haven't found it up to now, give up and return None
  return None

def load_ligand_pdb(ligprefix,folders) :

  pdbfile = locate_file(ligprefix+".pdb",folders)
  
  # Cannot do anything without a pdb-file so raise an exception
  if pdbfile is None :
    raise SetupError("Ligand file %s.pdb could not be found"%ligprefix)

  return pdbfile,simulationobjects.PDBFile(filename=pdbfile)

def prep_ligand(files,charge,ligobj12,folders,settings) :
  """ Prepare a ligand completely so that a ProtoMS
      templatefile and a box of water around the ligand
      exist
  
  Parameters
  ----------   
  ligprefix - the filename of the ligand, except the file extension
  charge - the net charge of the ligand
  files  - a dictionary of filenames and a simulationobjects.pdb object
  ligobj12 - a simulationobjects.pdb object consisting of the first two ligands merged
  folders - a set of folders where to look for the ligand files
  settings - additional settings, an ArgParse objecy

  Returns
  -------
  None
  """

  ligprefix = get_prefix(files["pdb"])

  # Try to locate all necessary files for a ligand
  files["prepi"] = locate_file(ligprefix+".prepi",folders)
  files["frcmod"] = locate_file(ligprefix+".frcmod",folders)
  files["zmat"] = locate_file(ligprefix+".zmat",folders)
  files["tem"] = locate_file(ligprefix+".tem",folders)
  files["wat"] = locate_file(ligprefix+"_box.pdb",folders)

  print "\nSetting up ligand: %s..."%files["pdb"]

  # Check to see if we have a template file
  if files["tem"] is None : 
    resnam = files["obj"].residues[1].name
    if files["prepi"] is None :
      # Here we need to run Antechamber
      print "Running antechamber. Please check the output carefully"
      files["prepi"] = tools.run_antechamber(files["pdb"],charge,resnam)
      print "Created prepi-file: %s"%files["prepi"]
    if files["frcmod"] is None :
      # Here we need to run parmchk
      print "Running parmchk. Please check the output carefully"
      files["frcmod"] = tools.run_parmchk(files["pdb"])
      print "Created frcmod-file: %s"%files["frcmod"]
     
    # By this stage we should have all necessary files to make the template file
    # Now we need to call a routine to generate the template file
    #files["obj"] = simulationobjects.pdb(filename=files["pdb"])
    files["tem"] = ligprefix+".tem"
    tem = tools.build_template(temfile=files["tem"],prepifile=files["prepi"],zmatfile=files["zmat"],
                        frcmodfile=files["frcmod"],resname=resnam)
    tem.write(files["tem"])
    if files["zmat"] is None :
      files["zmat"] = ligprefix+".zmat"
      tem.templates[0].write_zmat(files["zmat"])
      print "Created zmatrix (%s) for ligand. Please check the output carefully"%files["zmat"]
    print "Created ProtoMS template-file (%s) for ligand. Please check the output carefully"%files["tem"]
                          
  # Check to see if we have solvated the ligand
  if files["wat"] is None :
    # Here we need to call the solvate routine

    # Try to find a default water box
    if settings.waterbox is None :
      waterbox = simulationobjects.standard_filename("wbox_"+settings.watmodel+".pdb","data")
    else :
      waterbox = args.waterbox
    if not os.path.isfile(waterbox) : 
      raise SetupError("Could not find file (%s) with pre-equilibrated waters"%waterbox)

    # Setting the solute, either a pdb filename or a simulationobjects.pdb object
    if ligobj12 is None :
      solute = files["obj"]
    else :
      solute = ligobj12
    # Calling the routine
    print "Created waterbox-file: %s"%(ligprefix+"_box.pdb")
    files["wat"] = ligprefix+"_box.pdb"
    boxpdb = tools.solvate(waterbox, ligand=solute, protein=None,
                           geometry="box",padding=10.0, radius=30.0, center="cent",
                           namescheme="ProtoMS")
    boxpdb.write(files["wat"])

  return files
  
def prep_protein(protprefix,ligands,watprefix,folders,settings) :
  """ Prepare a protein completely such that a scoop pdb-file
      and a droplet of water exists
  
  Parameters
  ----------   
  protprefix - the filename of the protein, except the file extension
  ligands - a simulationobjects.pdb object with all ligands merged
  watprefix - the filename of the water sphere, except the file extension
  folders - a set of folders where to look for the ligand files
  settings - additional settings, an ArgParse objecy
  limit - minimum difference between number of residues in protein and scoop

  Returns
  -------
  the filenames of the protein and water droplet pdb-file
  """

  protprefix = get_prefix(protprefix)
  watprefix  = get_prefix(watprefix)

  # Try to locate necessary protein files
  protein_orig_file = locate_file(protprefix+".pdb",folders)
  protein_pms_file = locate_file(protprefix+"_pms.pdb",folders)
  protein_scoop_file = locate_file(protprefix+"_scoop.pdb",folders)
  protein_water = locate_file(watprefix+".pdb",folders)
  
  # Cannot do anything without a pdb-file, so raise an exception
  if protein_orig_file is None and protein_scoop_file is None and protein_pms_file is None:
    raise SetupError("Protein file (%s.pdb) and protein scoop file (%s_scoop.pdb) and protein pms file (%s_pms.pdb) could not be found"%(protprefix,protprefix,protprefix))

  print "\nSetting up protein: %s..."%protein_orig_file

  protobj = None
  if protein_scoop_file is None and protein_pms_file is None :
    # Start with an object for the original pdb-file
    protobj = simulationobjects.PDBFile(filename=protein_orig_file)
  
    # Trying to locate a default conversions file
    if settings.atomnames is None :
      conversionfile = simulationobjects.standard_filename("atomnamesmap.dat","data")
    else :
      conversionfile = args.atomnames
    if not os.path.isfile(conversionfile) : 
      raise SetupError("Could not find file (%s) with atom name conversions"%conversionfile)

    # Converting to ProtoMS atom names
    protobj = tools.pdb2pms(protobj,"amber",conversionfile)

    # Converting water molecules to specified model
    protobj = tools.convertwater(protobj,settings.watmodel)

    # Defining the center of the scoop...
    if ligands is None :
      if settings.center != None :
        ligobj = settings.center
      else :
        protobj.getCenter()
        ligobj = "%f %f %f" % tuple(protobj.center)
        print "Warning: No specified center for protein scoop. Using the center of the protein." 
    else :
      ligobj = ligands

    # Here we need to call the routine to make a scoop
    protobj_scooped = tools.scoop(protobj,ligobj,
                                  innercut=settings.innercut,outercut=settings.outercut,
                                  flexin=settings.flexin,flexout=settings.flexout)

    nresdiff = len(protobj_scooped.residues)-len(protobj.residues)
    print "Created scoop-pdb file by removing %d residues: %s"%(nresdiff,protein_scoop_file)
    if nresdiff < settings.scooplimit:
      print "Discarding scoop. Number of residues removed from the protein is too small (%d). Created %s_pms.pdb instead."%(nresdiff,protprefix)
      protein_pms_file = get_prefix(protein_orig_file)+"_pms.pdb"
      protobj.write(filename=protein_pms_file,header='REMARK Original file %s\nREMARK Atoms renamed according to ProtoMS naming standards.\n'%protein_orig_file)
    else :
      protobj = protobj_scooped
      protein_scoop_file = get_prefix(protein_orig_file)+"_scoop.pdb"
      protobj.write(protein_scoop_file, renumber=True)
      
  if protein_water is None :
    # Here we need to call the routine to solvate the protein

    # Try to find a default water box
    if settings.waterbox is None :
      waterbox = simulationobjects.standard_filename("wbox_"+settings.watmodel+".pdb","data")
    else :
      waterbox = args.waterbox
    if not os.path.isfile(waterbox) : 
      raise SetupError("Could not find file (%s) with pre-equilibrated waters"%waterbox)

    # Setting the solute, either a pdb filename or a simulationobjects.pdb object
    if protobj is None :
      solute = protein_scoop_file
    else :
      solute = protobj

    # Calling the routine
    protein_water = watprefix+".pdb"
    print "Created water cap-file: %s"%(protein_water)
    cappdb = tools.solvate(waterbox, ligand=ligands, protein=solute,
                          geometry="droplet",padding=10.0, radius=args.capradius, center="cent",
                          namescheme="ProtoMS")
    cappdb.write(protein_water)

  if protein_pms_file != None:
    return protein_pms_file,protein_water
  else:
    return protein_scoop_file,protein_water
  

if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program setup and run a ProtoMS simulations")
  parser.add_argument('-s','--simulation',choices=["none","equilibration","sampling","dualtopology"],help="the kind of simulation to setup",default="none")
  parser.add_argument('--dovacuum',action='store_true',help="turn on vacuum simulation for simulation types equilibration and sampling",default=False)
  parser.add_argument('-f','--folders',nargs="+",help="folders to search for files ",default=["."])
  parser.add_argument('-p','--protein',help="the prefix of the protein")
  parser.add_argument('-l','--ligand',nargs="+",help="the prefix of the ligand(s)")
  parser.add_argument('-w','--water',help="the prefix of the water/solvent",default="water")
  parser.add_argument('-c','--cmdfile',help="the prefix of the command file",default="run")
  parser.add_argument('--charge',nargs="+",type=float,help="the net charge of each ligand")
  parser.add_argument('--lambdas',nargs="+",type=float,help="the lambda values or the number of lambdas",default=[16])
  parser.add_argument('--center',help="the center of the scoop, if ligand is not available, either a string or a file with the coordinates",default=None)
  parser.add_argument('--innercut',type=float,help="maximum distance from ligand defining inner region of the scoop",default=16.0)
  parser.add_argument('--outercut',type=float,help="maximum distance from ligand defining outer region of the scoop",default=20.0)
  parser.add_argument('--flexin',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the inner region",default="flexible")
  parser.add_argument('--flexout',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the outer region",default="sidechain")
  parser.add_argument('--capradius',type=float,help="the radius of the droplet around the protein",default=30.0)
  parser.add_argument('--nequil',type=int,help="the number of equilibration steps",default=5E6)
  parser.add_argument('--nprod',type=int,help="the number of production steps",default=40E6)
  parser.add_argument('--dumpfreq',type=int,help="the output dump frequency",default=1E5)
  parser.add_argument('--atomnames',help="a file with atom name conversions")
  parser.add_argument('--waterbox',help="a file with pre-equilibrated water molecules")
  parser.add_argument('--watmodel',help="the name of the water model. Default = tip4p",choices=[ 'tip3p', 'tip4p'],default='tip4p')
  parser.add_argument('--scooplimit',help="the minimum difference between number of residues in protein and scoop for scoop to be retained",default=10)
  args = parser.parse_args()
  
  # Adds current folder to the folders
  args.folders.append(".")
  
  if args.protein is None and args.ligand is None :
    print "Nothing to do, so exit!"
    quit()
 
  # Set $PROTOMSHOME 
  if os.getenv("PROTOMSHOME") is None :
    str = os.path.dirname(os.path.abspath(__file__))
    print "Setting PROTOMSHOME to %s"%str
    os.environ["PROTOMSHOME"] = str # This does not change the original shell
     
  # Prepare each given ligand
  ligand_files = {} # This will be filled with a dictionary of filenames for each ligand
  ligands = [] # This will be a list of ligands
  ligpdbs = None # This will be a list of ligand pdb-files
  ligtems = None # This will be a list of ligand template files
  ligobjs = None # This will hold a merged pdb object of all ligand pdb objects
  ligand_water = None
  if args.ligand is not None :
    # Read in each ligand pdb file and create a pdb object
    for l in args.ligand :
      prefix = get_prefix(l)
      ligands.append(prefix)
      ligand_files[prefix] = {}
      ligand_files[prefix]["pdb"],ligand_files[prefix]["obj"] = load_ligand_pdb(prefix,args.folders)

    # Create merge pdb objects
    if len(ligands) >= 2 :
      ligobj12 = simulationobjects.merge_pdbs(ligand_files[l]["obj"] for l in ligands[:2])  
    else :
      ligobj12 = ligand_files[ligands[0]]["obj"]
    if len(ligands) > 1 :
      ligobjs = simulationobjects.merge_pdbs(ligand_files[l]["obj"] for l in ligands)
    else :
      ligobjs = ligand_files[ligands[0]]["obj"]
 
    # Now do the prepartions
    for i,l in enumerate(args.ligand) :
      if args.charge is not None and i < len(args.charge): 
        charge = args.charge[i]
      else :
        charge = 0
      if i > 1 : ligobj12 = None
      prefix = get_prefix(l)
      prep_ligand(ligand_files[prefix],charge,ligobj12,args.folders,args)

    ligpdbs = [ligand_files[l]["pdb"] for l in ligands]
    ligtems = [ligand_files[l]["tem"] for l in ligands]
    ligand_water = ligand_files[ligands[0]]["wat"]

    # Here we will merge ligand template files if there is more than one
    if len(ligtems) > 1 :
      temfile = tools.simulationobjects.TemplateFile(ligand_files[ligands[0]]["tem"])
      for l in ligands[1:] :
        temfile2 = tools.simulationobjects.TemplateFile(ligand_files[l]["tem"])
        temfile.append(temfile2)
      allnames = "-".join(t.name.lower() for t in temfile.templates)
      ligtems = [allnames+".tem"]
      temfile.write(allnames+".tem")
      print "\nCreated a re-numbered template files for all ligands: %s"%ligtems[0]
    
  # Prepare the protein
  protein_file = None
  water_file = None
  if args.protein is not None :
    protein_file,water_file = prep_protein(args.protein,ligobjs,args.water,args.folders,args)

  # Create ProtoMS command files
  free_cmd,bnd_cmd = tools.generate_input(protein_file,ligpdbs,ligtems,water_file,ligand_water,args)
  if free_cmd is not None : 
    free_cmd.writeCommandFile(args.cmdfile+"_free.cmd")
  if bnd_cmd is not None : 
    bnd_cmd.writeCommandFile(args.cmdfile+"_bnd.cmd")    
      
    
      
    


