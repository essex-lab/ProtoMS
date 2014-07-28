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
  
def _locate_file(filename,folders) :
  """ 
  Find a file 
  
  Tries to find the file as it is given or in
  any of the folders given

  Parameters
  ----------
  filename : string
    the name of the file to find
  folders : list of strings
    folders to search for the file

  Returns
  -------
  string or None
    the full filename or None if the file could not be found
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

def _merge_templates(templates) :
  """
  Merge template files
  
  Parameters
  ----------
  templates : list of string
    names of all the template files

  Returns
  -------
  list of strings
    modify list of names
  """
  temfile = tools.merge_templates(templates)
  allnames = "-".join(t.name.lower() for t in temfile.templates)
  templates = [allnames+".tem"]
  temfile.write(allnames+".tem")
  print "Created a re-numbered template files for all ligands: %s"%templates[0]
  return templates

def _load_ligand_pdb(ligprefix,folders) :
  """
  Load a ligand pdb file
  
  Parameters
  ----------
  ligprefix : string
    the filename of ligand pdb-file, without the extension
  folders : list of string
    folders to look for the pdb-file in
 
  Returns
  -------
  string
    the full filename of the pdb-file
  PDBFile
    an instance that hold the pdb-file structure

  Raises
  ------
  SetupError
    if the pdb-file cannot be found
  """
  pdbfile = _locate_file(ligprefix+".pdb",folders)
  
  # Cannot do anything without a pdb-file so raise an exception
  if pdbfile is None :
    raise simulationobjects.SetupError("Ligand file %s.pdb could not be found"%ligprefix)

  return pdbfile,simulationobjects.PDBFile(filename=pdbfile)

def _prep_ligand(files,charge,ligobj12,folders,settings) :
  """ 
  Prepare a ligand completely so that a ProtoMS
  templatefile and a box of water around the ligand exist
  
  Parameters
  ----------   
  files : dictionary
    filenames and a PDBFile object associated with this ligand
  charge : int
    the net charge of the ligand
  ligobj12 : PDBFile
    an instance consisting of the first two ligands merged
  folders : list of string
    folders where to look for the ligand files
  settings : Namespace (from argparse) 
    additional settings

  Returns
  -------
  dictionary
    filenames and a PDBFile object associates with this ligand

  Raises
  ------
  SetupError
    if a default water box could not be found
  """

  ligprefix = _get_prefix(files["pdb"])

  # Try to locate all necessary files for a ligand
  files["prepi"] = _locate_file(ligprefix+".prepi",folders)
  files["frcmod"] = _locate_file(ligprefix+".frcmod",folders)
  files["zmat"] = _locate_file(ligprefix+".zmat",folders)
  files["tem"] = _locate_file(ligprefix+".tem",folders)
  files["wat"] = _locate_file(ligprefix+"_box.pdb",folders)

  print "\nSetting up ligand: %s..."%files["pdb"]

  # Check to see if we have a template file
  if files["tem"] is None : 
    resnam = files["obj"].residues[1].name # Set the prepi name and template name to the residue name
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

    # Setting the solute
    if ligobj12 is None :
      solute = files["obj"]
    else :
      solute = ligobj12
    # Calling the routine
    files["wat"] = ligprefix+"_box.pdb"
    print "Created waterbox-file: %s"%(files["wat"])
    boxpdb = tools.solvate(settings.waterbox, ligand=solute, protein=None,
                           geometry="box",padding=10.0, radius=30.0, center="cent",
                           namescheme="ProtoMS")
    boxpdb.write(files["wat"])

  return files
  
def _prep_protein(protprefix,ligands,watprefix,folders,settings) :
  """ 
  Prepare a protein completely such that a scoop pdb-file
  and a droplet of water exists
  
  Parameters
  ----------   
  protprefix : string 
    the filename of the protein
  ligands : PDBFile
    an instance with all ligands merged
  watprefix : string 
    the filename of the water sphere
  folders : list of strings
    folders where to look for the protein files
  settings : Namespace (from argparse) 
    additional settings

  Returns
  -------
  string
    the filename of the protein pdb-file to use
  string
    the filename of the water droplet

  Raises
  ------
  SetupError
    if protein pdb-file cannot be found
    if default atom name conversion file cannot be found
    if default water box cannot be found
  """

  protprefix = _get_prefix(protprefix)
  watprefix  = _get_prefix(watprefix)

  # Try to locate necessary protein files
  protein_orig_file = _locate_file(protprefix+".pdb",folders)
  protein_pms_file = _locate_file(protprefix+"_pms.pdb",folders)
  protein_scoop_file = _locate_file(protprefix+"_scoop.pdb",folders)
  protein_water = _locate_file(watprefix+".pdb",folders)
  
  # Cannot do anything without a pdb-file, so raise an exception
  if protein_orig_file is None and protein_scoop_file is None and protein_pms_file is None:
    raise simulationobjects.SetupError("Protein file (%s.pdb) and protein scoop file (%s_scoop.pdb) and protein pms file (%s_pms.pdb) could not be found"%(protprefix,protprefix,protprefix))

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
      raise simulationobjects.SetupError("Could not find file (%s) with atom name conversions"%conversionfile)

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

    nresdiff = len(protobj.residues)-len(protobj_scooped.residues)
    protein_scoop_file = _get_prefix(protein_orig_file)+"_scoop.pdb"
    print "Created scoop-pdb file by removing %d residues: %s"%(nresdiff,protein_scoop_file)
    if nresdiff < settings.scooplimit:
      protein_pms_file = _get_prefix(protein_orig_file)+"_pms.pdb"
      print "Discarding scoop. Number of residues removed from the protein is too small (%d). Created %s instead."%(nresdiff,protein_pms_file)
      protobj.write(filename=protein_pms_file,header='REMARK Original file %s\nREMARK Atoms renamed according to ProtoMS naming standards.\n'%protein_orig_file)
      protein_scoop_file = None
    else :
      protobj = protobj_scooped
      protobj.write(protein_scoop_file, renumber=True)
      
  if protein_water is None :
   
    if protobj is None :
      if protein_scoop_file is not None :
        solute = protein_scoop_file
      else :
        solute = protein_pms_file
    else :
      solute = protobj

    # Calling the routine
    protein_water = watprefix+".pdb"
    print "Created water cap-file: %s"%(protein_water)
    cappdb = tools.solvate(settings.waterbox, ligand=ligands, protein=solute,
                          geometry="droplet",padding=10.0, radius=args.capradius, center="cent",
                          namescheme="ProtoMS")
    cappdb.write(protein_water)

  if protein_pms_file is not None:
    return protein_pms_file,protein_water
  else:
    return protein_scoop_file,protein_water

def _prep_singletopology(pdbs,templates1,settings) :
  """
  Prepare templates for single topology
  
  Parameters
  ----------
  pdbs : list of strings
    names of ligand pdb files
  templares1 : list of strings
    names of ligand template files
  settings : Namespace (from argparse) 
    additional settings
    
  Returns
  -------
  list of string
    the template files for electrostatic leg
  list of string
    the template files for the van der Waals leg
  """
  tem1 = simulationobjects.TemplateFile(templates1[0])
  tem2 = simulationobjects.TemplateFile(templates1[1])

  print "\nSetting up single-topology correspondance map and templates..."
  eletem,vdwtem,cmap = tools.make_single(tem1,tem2,pdbs[0],pdbs[1],settings.singlemap)
  
  pdbs.pop(1)
  templates1.pop(1)
  templates2 = list(templates1)

  prefix = os.path.join(os.path.dirname(pdbs[0]),tem1.templates[0].name.lower()+"t"+tem2.templates[0].name.lower())
  templates1[0] = prefix+"_ele.tem" 
  templates2[0] = prefix+"_vdw.tem"

  eletem.write(templates1[0])
  vdwtem.write(templates2[0])
  print "\nCreated template %s for electrostatic perturbation. Please check the output carefully."%templates1[0]
  print "Created template %s for van der Waals-perturbation. Please check the output carefully."%templates2[0]
  
  if args.singlemap is None : settings.singlemap = "single_cmap.dat"
  print "Saved correspondance map to: %s"%settings.singlemap
  tools.write_map(cmap,settings.singlemap)

  return templates1,templates2

def _prep_gcmc(ligands,ligand_files,settings) :
  """
  Prepare water box for GCMC calculations
  
  Parameters
  ----------
  ligands : list of strings
    names of all ligands loaded
  ligand_files : dictionary
    files and objects associated with each ligand
  settings : Namespace (from argparse) 
    additional settings
    
  Returns
  -------
  string
    the filename of the created water box
  """
  
  def pdb2box(pdbobj) :
    BOX_PADDING = 2.0
    # Create a box around the solute and pad it with two Angstromgs
    box = pdbobj.getBox()
    box["origin"] = box["origin"] - BOX_PADDING
    box["len"] = box["len"] + 2.0*BOX_PADDING
    # Save it to disc
    boxpdb = "gcmc_box.pdb"
    simulationobjects.write_box(boxpdb,box)
    print "\nCreated %s to visualize GCMC simulation box. Please check the output carefully"%boxpdb
    return boxpdb
    
  if ligands :
    # Use the ligand to define the GCMC box
    boxpdb = pdb2box(ligand_files[ligands[0]]["obj"])  
  else :
    if settings.gcmcbox is not None :
      gcmcboxobj = simulationobjects.PDBFile(filename=settings.gcmcbox) 
      # Try the find the box in the header of the file
      box = simulationobjects.find_box(gcmcboxobj)
      if box is None :
        # Else take it as a ligand
       boxpdb = pdb2box(gcmcboxobj)
      else :
        boxpdb = settings.gcmcbox
    else :
      simulationobjects.SetupError("Cannot define a GCMC simulation box without a ligand and without the gcmcbox setting") 
    
  # Use the flood option in solvate
  boxobj = tools.solvate(settings.waterbox, ligand=boxpdb, protein=None,
                         geometry="flood",namescheme="ProtoMS")
  # Need to change the residue name of the created waters
  for sol in boxobj.solvents :
    for atom in boxobj.solvents[sol].atoms : atom.resname = "WAT"
  # Write the box to disc
  boxname = "gcmc_wat.pdb"
  boxobj.write(boxname) 
  print "\nCreated %s; it contains the GCMC simulation waters. Please check the output carefully"%boxname       
  return boxname

  
if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program setup and run a ProtoMS simulations")
  parser.add_argument('-s','--simulation',choices=["none","equilibration","sampling","dualtopology","singletopology","gcmc"],help="the kind of simulation to setup",default="none")
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
  parser.add_argument('--adams',nargs="+",type=float,help="the Adam/B values for the GCMC",default=0)
  parser.add_argument('--gcmcwater',help="a pdb file with a box of water to do GCMC on")
  parser.add_argument('--gcmcbox',help="a pdb file with box dimensions for the GCMC box")
  parser.add_argument('--singlemap',help="the correspondance map for single-topology")
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

  # Try to find a default water box
  if args.waterbox is None :
    args.waterbox = simulationobjects.standard_filename("wbox_"+args.watmodel.lower()+".pdb","data")
  if not os.path.isfile(args.waterbox) : 
    raise simulationobjects.SetupError("Could not find file (%s) with pre-equilibrated waters"%waterbox)
     
  # Prepare each given ligand
  ligand_files = {} # This will be filled with a dictionary of filenames for each ligand
  ligands = [] # This will be a list of ligands
  ligpdbs = None # This will be a list of ligand pdb-files
  ligtems = None # This will be a list of ligand template files
  ligobjs = None # This will hold a merged pdb object of all ligand pdb objects
  ligand_water = None # This will hold the filename of the free-leg waterbox
  if args.ligand is not None :
    # Read in each ligand pdb file and create a pdb object
    for l in args.ligand :
      prefix = _get_prefix(l)
      ligands.append(prefix)
      ligand_files[prefix] = {}
      ligand_files[prefix]["pdb"],ligand_files[prefix]["obj"] = _load_ligand_pdb(prefix,args.folders)

    # Create merge pdb objects
    if len(ligands) >= 2 :
      ligobj12 = simulationobjects.merge_pdbs(ligand_files[l]["obj"] for l in ligands[:2])  
    else :
      ligobj12 = ligand_files[ligands[0]]["obj"]
    if len(ligands) > 1 :
      ligobjs = simulationobjects.merge_pdbs(ligand_files[l]["obj"] for l in ligands)
    else :
      ligobjs = ligand_files[ligands[0]]["obj"]
 
    # Now do the preparations
    for i,l in enumerate(args.ligand) :
      if args.charge is not None and i < len(args.charge): 
        charge = args.charge[i]
      else :
        charge = 0
      if i > 1 : ligobj12 = None
      prefix = _get_prefix(l)
      _prep_ligand(ligand_files[prefix],charge,ligobj12,args.folders,args)

    ligpdbs = [ligand_files[l]["pdb"] for l in ligands]
    ligtems = [ligand_files[l]["tem"] for l in ligands]
    ligand_water = ligand_files[ligands[0]]["wat"]
 
    # Here we need to make single topology templates, if requested
    if args.simulation == "singletopology" :
      ligtems,ligtems2 = _prep_singletopology(ligpdbs,ligtems,args)

    # Here we will merge ligand template files if there is more than one
    if len(ligtems) > 1 :
      print ""
      ligtems = _merge_templates(ligtems)
      if args.simulation == "singletopology" : ligtems2 = _merge_templates(ligtems2)    
    
  # Prepare the protein
  protein_file = None
  water_file = None
  if args.protein is not None :
    protein_file,water_file = _prep_protein(args.protein,ligobjs,args.water,args.folders,args)

  # Extra preparation for GCMC
  if args.simulation == "gcmc" and args.gcmcwater is None:
    args.gcmcwater = _prep_gcmc(ligands,ligand_files,args)    
      
  # Create ProtoMS command files
  postfix = ""
  if args.simulation == "singletopology" : 
    postfix = "_ele"
    setattr(args,"outfolder","out_ele")
  free_cmd,bnd_cmd = tools.generate_input(protein_file,ligpdbs,ligtems,water_file,ligand_water,args)
  if free_cmd is not None : 
    free_cmd.writeCommandFile(args.cmdfile+postfix+"_free.cmd")
  if bnd_cmd is not None : 
    bnd_cmd.writeCommandFile(args.cmdfile+postfix++"_bnd.cmd")    
  
  if args.simulation == "singletopology" :
    setattr(args,"outfolder","out_vdw")
    free_cmd,bnd_cmd = tools.generate_input(protein_file,ligpdbs,ligtems2,water_file,ligand_water,args)
    if free_cmd is not None : 
      free_cmd.writeCommandFile(args.cmdfile+"_vdw_free.cmd")
    if bnd_cmd is not None : 
      bnd_cmd.writeCommandFile(args.cmdfile+"_vdw_bnd.cmd")      
      
    
      
    


