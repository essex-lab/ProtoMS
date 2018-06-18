#!/usr/bin/env python2.7
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
import sys
import subprocess
import logging
import time

import numpy as np

import tools
from tools import simulationobjects

# Logger is used globally
logger = simulationobjects.setup_logger("protoms_py.log")

def _is_float(num) :
  """
  Check whether a string is convertible to a float
  
  Parameters
  ----------
  num : string
    the string which might be convertible to float
  
  Returns
  -------
  boolean
    whether the string is convertible to float
  """
  try :
    float(num)
  except :
    return False
  return True

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

def _merge_templates(templates,tarlist) :
  """
  Merge template files
  
  Parameters
  ----------
  templates : list of string
    names of all the template files
  tarlist : list of string
    name of files that can be stored away

  Returns
  -------
  list of strings
    modify list of names
  """
  for tem in templates : tarlist.append(tem) # All of the original templates can safely be stored away
  temfile = tools.merge_templates(templates)
  if isinstance(temfile,simulationobjects.TemplateFile):
    allnames = "-".join(t.name.lower() for t in temfile.templates)
    templates = [allnames+".tem"]
    temfile.write(allnames+".tem")
    logger.info("Created a re-numbered template files for all ligands: %s"%templates[0])
    return templates
  else:
    return [temfile]

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
    msg = "Ligand file %s.pdb could not be found"%ligprefix
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  return pdbfile,simulationobjects.PDBFile(filename=pdbfile)

def _prep_ligand(files,first,charge,ligobj12,folders,tarlist,settings) :
  """ 
  Prepare a ligand completely so that a ProtoMS
  templatefile and a box of water around the ligand exist
  
  Parameters
  ----------   
  files : dictionary
    filenames and a PDBFile object associated with this ligand
  first : bool
    flag to indicate if the ligand is the first one
  charge : int
    the net charge of the ligand
  ligobj12 : PDBFile
    an instance consisting of the first two ligands merged
  folders : list of string
    folders where to look for the ligand files
  tarlist : list of string
    name of files that can be stored away
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
  files["wat"] = _locate_file(ligprefix+"_box.pdb",folders)

  logger.info("")
  logger.info("Setting up ligand: %s..."%files["pdb"])

  if len(files["obj"].residues) < 1 :
    if len(files["obj"].solvents) > 0 :
      raise simulationobjects.SetupError("Your ligand in %s is recognized as solvent. Please change the residue name." %files["obj"].name)
    else :
      raise simulationobjects.SetupError("No residues found in %s."%files["obj"].name)

  # Get the ligand name from the pdb header
  if 'HEADER' in files["obj"].header:
    words = files["obj"].header.strip().split()
    ligname = words[words.index('HEADER')+1]
  else:
    logger.info('Unable to find header information in the PDB-file, adding it automatically.')
    ligname = files["obj"].residues[1].name
    files["obj"].header = files["obj"].header + "HEADER " + ligname + "\n"
    files["obj"].write(files["pdb"])

  # Try to locate template for the ligand
  files["tem"] = None
  if settings.template is None:
    files["tem"] = _locate_file(ligprefix+".tem",folders)
  else:
    for f in settings.template :
      tempfile = _locate_file(f,folders)
      if tempfile is None : 
        files["tem"] = tempfile
      else :
        tem = simulationobjects.TemplateFile(filename=tempfile)
        if ligname.lower() in [t.name.lower() for t in tem.templates]:
          files["tem"] = tempfile
          break

  # Check to see if we have a template file
  if files["tem"] is None : 
    resnam = files["obj"].residues[1].name # Set the prepi name and template name to the residue name
    if files["prepi"] is None :
      # Here we need to run Antechamber
      logger.info("Running antechamber. Please check the output carefully")
      files["prepi"] = tools.run_antechamber(files["pdb"],charge,resnam)
      logger.info("Created prepi-file: %s"%files["prepi"]) 
      tarlist.append(files["prepi"])
    if files["frcmod"] is None :
      # Here we need to run parmchk
      logger.info("Running parmchk. Please check the output carefully")
      files["frcmod"] = tools.run_parmchk(files["pdb"])
      logger.info("Created frcmod-file: %s"%files["frcmod"])
      tarlist.append(files["frcmod"])
     
    # By this stage we should have all necessary files to make the template file
    files["tem"] = ligprefix+".tem"
    tem = tools.build_template(temfile=files["tem"],prepifile=files["prepi"],zmatfile=files["zmat"],
                        frcmodfile=files["frcmod"],resname=resnam,gaffversion=settings.gaff)
    tem.write(files["tem"])
    if files["zmat"] is None :
      files["zmat"] = ligprefix+".zmat"
      tem.templates[0].write_zmat(files["zmat"])
      logger.info("Created zmatrix (%s) for ligand. Please check the output carefully"%files["zmat"])
      tarlist.append(files["zmat"])
    logger.info("Created ProtoMS template-file (%s) for ligand. Please check the output carefully"%files["tem"])
    
  # Check to see if we have solvated the ligand
  if files["wat"] is None :

    # Setting the solute
    if ligobj12 is None :
      solute = files["obj"]
    else :
      solute = ligobj12
    # Calling the routine
    files["wat"] = ligprefix+"_box.pdb"
    boxpdb = tools.solvate(settings.waterbox, ligand=solute, protein=None,
                           geometry="box",padding=10.0, radius=30.0, center="cent",
                           namescheme="ProtoMS")
    
    boxpdb.write(files["wat"])
    logger.info("Created waterbox-file: %s"%(files["wat"]))
    if (settings.simulation in ["dualtopology","singletopology"] and not first) or (settings.simulation in ["gcmc","jaws1","jaws2"]) :
      tarlist.append(files["wat"])

  return files
  
def _prep_protein(protprefix,ligands,watprefix,folders,tarlist,settings) :
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
  tarlist : list of string
    name of files that can be stored away
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

  if protprefix is not None :
    protprefix = _get_prefix(protprefix)
  else :
    protprefix = _get_prefix(settings.scoop)
  watprefix  = _get_prefix(watprefix)
  if settings.scoop is None :
    scoopprefix = protprefix+"_scoop"
  else :
    scoopprefix = _get_prefix(settings.scoop)

  # Try to locate necessary protein files
  protein_orig_file = _locate_file(protprefix+".pdb",folders)
  protein_pms_file = _locate_file(protprefix+"_pms.pdb",folders)
  protein_scoop_file = _locate_file(scoopprefix+".pdb",folders)
  protein_water = _locate_file(watprefix+".pdb",folders)
  
  # Cannot do anything without a pdb-file, so raise an exception
  if settings.protein is None and protein_scoop_file is None:
    msg = "Specified scoop (%s.pdb) could not be found"%scoopprefix
    logger.error(msg)
    raise simulationobjects.SetupError(msg)
  if protein_orig_file is None and protein_scoop_file is None and protein_pms_file is None:
    msg = "Protein file (%s.pdb) and protein scoop file (%s.pdb) and protein pms file (%s_pms.pdb) could not be found"%(protprefix,scoopprefix,protprefix)
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  logger.info("")
  logger.info("Setting up protein: %s..."%protein_orig_file)

  if settings.scoop is not None and protein_scoop_file is None:
    logger.info("Specified scoop (%s.pdb) not found. Ignoring..."%scoopprefix)

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
      msg = "Could not find file (%s) with atom name conversions"%conversionfile
      logger.error(msg)
      raise simulationobjects.SetupError(msg)

    # Converting to ProtoMS atom names
    protobj2 = tools.pdb2pms(protobj,"amber",conversionfile)
    abortedconv = protobj2 == protobj  # Indicates abortion of conversion
    protobj = protobj2

    # Converting water molecules to specified model
    protobj = tools.convertwater(protobj,settings.watmodel)

    # Defining the center of the scoop...
    if ligands is None and settings.gcmcbox is None:
      if settings.center != None :
        ligobj = settings.center
      else :
        protobj.getCenter()
        ligobj = "%f %f %f" % tuple(protobj.center)
        logger.warning("Warning: No specified center for protein scoop. Using the center of the protein." )
    elif ligands is None and settings.gcmcbox is not None :
      ligobj = simulationobjects.PDBFile(filename=settings.gcmcbox[0]) 
      logger.info("Scooping protein around %s..." %  settings.gcmcbox[0])
    else :
      ligobj = ligands

    # Here we need to call the routine to make a scoop
    protobj_scooped = tools.scoop(protobj,ligobj,
                                  innercut=settings.innercut,outercut=settings.outercut,
                                  flexin=settings.flexin,flexout=settings.flexout,scooplimit=settings.scooplimit)

    nresdiff = len(protobj.residues)-len(protobj_scooped.residues)
    protein_scoop_file = _get_prefix(protein_orig_file)+"_scoop.pdb"
    #logger.info("Created scoop-pdb file by removing %d residues: %s"%(nresdiff,protein_scoop_file))  
    protobj = protobj_scooped
    if nresdiff < settings.scooplimit:
      if not abortedconv : 
        protein_pms_file = _get_prefix(protein_orig_file)+"_pms.pdb"
        logger.info("Created %s instead."%protein_pms_file)
        header = protobj.header
        header += 'REMARK Atoms renamed according to ProtoMS naming standards.\n'
        protobj.write(filename=protein_pms_file,header=header, solvents = False)
      else : 
        #logger.info("Not scooping. Number of residues removed from the protein is too small.")
        protein_pms_file = protein_orig_file
      tarlist.append(protein_scoop_file)
      protein_scoop_file = None
    else :
      logger.info("Created scoop-pdb file by removing %d residues: %s"%(nresdiff,protein_scoop_file))  
      protobj.write(protein_scoop_file, renumber=True, solvents = False)
      
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
    cappdb = tools.solvate(settings.waterbox, ligand=ligands, protein=solute,
                          geometry="droplet",padding=10.0, radius=args.capradius, center="cent",
                          namescheme="ProtoMS")
    cappdb.write(protein_water)
    logger.info("Created water cap-file: %s"%(protein_water))

  if protein_pms_file is not None:
    return protein_pms_file,protein_water
  else:
    return protein_scoop_file,protein_water

def _prep_singletopology(pdbs,templates1,tarlist,settings) :
  """
  Prepare templates for single topology
  
  Parameters
  ----------
  pdbs : list of strings
    names of ligand pdb files
  templates1 : list of strings
    names of ligand template files
  tarlist : list of string
    name of files that can be stored away
  settings : Namespace (from argparse) 
    additional settings
    
  Returns
  -------
  list of string
    the template files for electrostatic leg
  list of string
    the template files for the van der Waals leg
  """
  try:
    tem1 = simulationobjects.TemplateFile(templates1[0])
    tem2 = simulationobjects.TemplateFile(templates1[1])
  except IndexError:
    msg = "Single topology calculations require two ligand templates to perturb between.\nYou only have one: %s \nIf you are trying to run an absolute free energy calculation\n(i.e. perturb to nothing), please use dual topology." %templates1[0]
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  # The original templates file can now be stored away
  tarlist.append(templates1[0])
  tarlist.append(templates1[1])

  logger.info("")
  logger.info("Setting up single-topology correspondance map and templates...")
  eletem,vdwtem,combtem,cmap = tools.make_single(tem1,tem2,pdbs[0],pdbs[1],settings.singlemap,gaffversion=settings.gaff)
  tools.summarize(eletem,vdwtem,logger.debug)  

  pdbs.pop(1)
  templates1.pop(1)
  templates2 = list(templates1)
  templates3 = list(templates1)

  prefix = os.path.join(os.path.dirname(pdbs[0]),tem1.templates[0].name.lower()+"t"+tem2.templates[0].name.lower())
  templates1[0] = prefix+"_ele.tem" 
  templates2[0] = prefix+"_vdw.tem"
  templates3[0] = prefix+"_comb.tem"

  eletem.write(templates1[0])
  vdwtem.write(templates2[0])
  combtem.write(templates3[0])
  logger.info("")
  logger.info("Created template %s for electrostatic perturbation. Please check the output carefully."%templates1[0])
  logger.info("Created template %s for van der Waals-perturbation. Please check the output carefully."%templates2[0])
  logger.info("Created template %s for combined perturbation. Please check the output carefully."%templates3[0])
  
  if args.singlemap is None : 
    settings.singlemap = "single_cmap.dat"
    tarlist.append(settings.singlemap)
  tools.write_map(cmap,settings.singlemap)
  logger.info("Saved correspondance map to: %s"%settings.singlemap)

  return templates1,templates2,templates3

def _prep_gcmc(ligands,ligand_files,waters,tarlist,settings) :
  """
  Prepare water box for GCMC or JAWS-1 calculations
  
  Parameters
  ----------
  ligands : list of strings
    names of all ligands loaded
  ligand_files : dictionary
    files and objects associated with each ligand
  waters : string
    the name of the filename containing solvation waters
  tarlist : list of string
    name of files that can be stored away
  settings : Namespace (from argparse) 
    additional settings
    
  Returns
  -------
  string
    the filename of the created water box
  """
  
  def pdb2box(pdbobj,padding=2.0) :
    boxpdb = "%s_box.pdb"%settings.simulation
    tools.make_gcmcbox(pdbobj,boxpdb,padding)
    logger.info("")
    logger.info("Created %s to visualize GCMC/JAWS-1 simulation box. Please check the output carefully"%boxpdb)
    return boxpdb

  def arrange_wats(box,water_file,out_name) :
    # find smallest box including all oxygens in my water_file
    waterobj = simulationobjects.PDBFile(filename = water_file)
    if len(waterobj.residues) !=0 :
      waterbox = waterobj.getBox(atomlist=[waterobj.residues.itervalues().next().atoms[0].name])
    else:
      waterbox = waterobj.getBox(atomlist=[waterobj.solvents.itervalues().next().atoms[0].name])

    # check whether the box of the water oxygens is within the box provided as argument
    box_origen_below = np.all(box['origin']< waterbox['origin'])
    box_end_avobe = np.all(box['origin']+box['len'] > waterbox['origin']+waterbox['len'])
    if not (box_origen_below and box_end_avobe) or (np.all(waterbox['len'] < 0.5) and len(waterobj.solvents) > 1) :
      # re-distribute the waters if required according to provided box
      logger.info("\nRedistributing your GCMC / JAWS1 waters according to the specified box\n")
      arranged_obj = tools.distribute_particles(box,water_file)
      arranged_obj.write(filename=out_name)
      return arranged_obj, out_name
    else :
      # make sure the header in the waterobj is the box dimensions which will be used in this simulation
      box_extremes = tuple(box['origin']) + tuple(np.add(box['origin'],box['len']))
      waterobj.header = "HEADER box %.3f %.3f %.3f %.3f %.3f %.3f \n" %box_extremes
      return waterobj, out_name

  def fill_box (settings,boxpdb,box,ghost_name) :
    # Fill the box with waters
    if settings.gcmcwater is None :
      # Creating ghost waters with a density 1.5 times the density of bulk water.
      box_volume = box['len'][0]*box['len'][1]*box['len'][2] 
      ghostnum = str(int(np.ceil(box_volume*0.0334*1.5)))     # 0.0334 waters per Angs.^3 is the number density of bulk water.
      ghostobj = tools.distribute_particles(box=box,particles=ghostnum,watermodel=settings.watmodel)
      ghostobj.write(filename=ghost_name)
#      ghostobj = tools.solvate(settings.waterbox, ligand=boxpdb, protein=None,
#                         geometry="flood",namescheme="ProtoMS")
#      boxinfo = tuple(box['origin']) + tuple(box['origin'] + box['len'])
#      ghostobj.header = "HEADER box %.4f %.4f %.4f %.4f %.4f %.4f\n" %boxinfo
#      for sol in ghostobj.solvents :
#        for atom in ghostobj.solvents[sol].atoms : atom.resname = "WA1"
    elif settings.gcmcwater.isdigit() :
      print settings.gcmcwater
      ghostobj = tools.distribute_particles(box,settings.gcmcwater,watermodel=settings.watmodel)
      ghostobj.write(filename=ghost_name)
    else :
      ghostobj, ghost_name = arrange_wats(box,settings.gcmcwater,ghost_name)
    return ghostobj, ghost_name    


  # Check consistency of command-line arguments
  if settings.gcmcwater is not None and not settings.gcmcwater.isdigit():
    gcmcwater = _locate_file(settings.gcmcwater,settings.folders)
    if gcmcwater is None :
      msg = "File %s given as gcmcwater could not be found."%settings.gcmcwater
      logger.error(msg)
      raise simulationobjects.SetupError(msg)
    settings.gcmcwater = gcmcwater
  if settings.gcmcbox is not None and len(settings.gcmcbox) is 1:
    gcmcbox = _locate_file(settings.gcmcbox[0],settings.folders)
    if gcmcbox is None :
      msg = "File %s given as gcmcbox could not be found."%settings.gcmcbox[0]
      logger.error(msg)
      raise simulationobjects.SetupError(msg)
    settings.gcmcbox = gcmcbox
  elif settings.gcmcbox is not None and len(settings.gcmcbox) < 6 :
    msg = "6 arguments expected to define the GCMC/JAWS1 box dimensions, %d provided: %s"%(len(settings.gcmcbox)," ".join(settings.gcmcbox))
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  ghost_name = "%s_wat.pdb"%settings.simulation
  write = True

  # If the gcmcbox has been set as a file
  if settings.gcmcbox is not None  and isinstance(settings.gcmcbox,str) :
    gcmcboxobj = simulationobjects.PDBFile(filename=settings.gcmcbox) 
    # Try to find the box in the header of the file
    box = simulationobjects.find_box(gcmcboxobj)
    if box is None :
      boxpdb = pdb2box(gcmcboxobj,padding=0.0)
      gcmcboxobj = simulationobjects.PDBFile(filename=boxpdb)
      box = simulationobjects.find_box(gcmcboxobj)
    else :
      boxpdb = settings.gcmcbox
    if "center" in box :
      box['origin'] = np.array([coord-box["len"][ind]/2 for ind,coord in enumerate(box["center"])])

    # Fill the box with waters
    ghostobj, ghost_name = fill_box (settings,boxpdb,box,ghost_name)
        

  # If the dimensions of the gcmc box have been provided
  elif settings.gcmcbox is not None and len(settings.gcmcbox) is 6 and [_is_float(value) for value in settings.gcmcbox] :
    # Generate the box pdb
    boxpdb = "%s_box.pdb"%settings.simulation
    try :
      box={}
      box['origin'] = np.array([float(value) for value in settings.gcmcbox[:3]])
      box['len'] = np.array([float(value) for value in settings.gcmcbox[3:]])
    except :
      raise simulationobjects.SetupError("The box dimensions %s could not be correctly interpreted"%settings.gcmcbox)
    simulationobjects.write_box(boxpdb,box)
    logger.info("")
    logger.info("Created %s to visualize GCMC/JAWS-1 simulation box. Please check the output carefully"%boxpdb)
    # Fill the box with waters
    ghostobj, ghost_name = fill_box (settings,boxpdb,box,ghost_name)
  
  # If no box information has been provided but we have a gcmcwater file
  elif settings.gcmcbox is None and settings.gcmcwater is not None and not settings.gcmcwater.isdigit() :
    waterobj = simulationobjects.PDBFile(filename = settings.gcmcwater)
    box = simulationobjects.find_box(waterobj)
    if box is None :
      raise simulationobjects.SetupError("Cannot set up simulation without a box.")
    else :
      ghostobj = waterobj
      ghost_name = waterobj.name
      write = False

  # If no information on the box has been provided
  elif ligands :
    # Use the ligand to define the GCMC box
    boxpdb = pdb2box(ligand_files[ligands[0]]["obj"])  
    gcmcboxobj = simulationobjects.PDBFile(filename=boxpdb)
    box = simulationobjects.find_box(gcmcboxobj)
    if "center" in box :
      box['origin'] = np.array([coord-box["len"][ind]/2 for ind,coord in enumerate(box["center"])])
    # Fill the box with waters
    ghostobj, ghost_name = fill_box (settings,boxpdb,box,ghost_name)
  else : 
    msg = "Cannot define a GCMC or JAWS-1 simulation box without a ligand and without the gcmcbox setting"
    logger.error(msg)
    simulationobjects.SetupError(msg) 
    
  # Write the waters to disc
  if write :
    ghostobj.write(ghost_name) 
    logger.info("")
    logger.info("Created %s; it contains the GCMC or JAWS-1 simulation waters. Please check the output carefully"%ghost_name)
  
  # Clear the GCMC/JAWS-1 box from solvation waters
  logger.info("")
  nrem,waters2 = tools.clear_gcmcbox(ghostobj,waters)
  if nrem > 0 :
    waters2_name = _get_prefix(str(waters2))+"_clr.pdb"
    logger.info("Created water cap-file: %s"%waters2_name)
    tarlist.append(waters)
    waters2.write(waters2_name)
    waters = waters2_name
  
  return ghost_name,waters

def _prep_jaws2(water_file,tarlist,settings) :
  """
  Prepare files for JAWS-2 simulation

  Parameters
  ----------
  water_file : string
    the name of the filename containing solvation waters
  tarlist : list of string
    name of files that can be stored away
  settings : Namespace (from argparse) 
    additional settings

  Returns
  -------
  PDBSet
    the set of single water molecule PDB files
  PDBSet
    the set of other water molecule PDB files
  string
    the name of the solvation waters

  Raises
  ------
  SetupError
    no gcmcwater specified in the settings
  """

  if settings.gcmcwater is None :
    msg = "You must set gcmcwater settings when preparing JAWS-2 input"
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  single_wat,other_wat = tools.split_waters(settings.gcmcwater)
  single_wat = tools.set_jaws2_box(single_wat)
  logger.info("")
  logger.info("Creating water PDB-files for JAWS-2 called jaws2_wat*.pdb and jaws2_not*.pdb")
  single_wat.write(["jaws2_wat%d.pdb"%(i+1) for i in range(len(single_wat.pdbs))])
  other_wat.write(["jaws2_not%d.pdb"%(i+1) for i in range(len(single_wat.pdbs))])

  nrem = 0
  for count,watobj in enumerate(single_wat.pdbs) :
    # Clear the JAWS-2 box from solvation waters
    watobj.name = "jaws2_wat%d.pdb"%(count+1)
    n,water_file = tools.clear_gcmcbox(watobj,water_file)
    nrem = nrem + n
  if nrem > 0 :
    waters2_name = _get_prefix(str(water_file))+"_clr.pdb"
    logger.info("Created water cap-file: %s"%waters2_name)
    tarlist.append(water_file)
    water_file.write(waters2_name)
    water_file = waters2_name
  else :
    water_file = water_file.name
 
  return single_wat,other_wat,water_file

def _cleanup(tarlist) :
  """
  Clean up extra files
  
  Parameters
  ----------
  tarlist : list of string
    the files to be cleaned up
  """
  tarlist2 = []
  for filename in tarlist :
    if filename in tarlist2 : continue
    if filename.find(os.environ["PROTOMSHOME"]) == 0 : continue
    tarlist2.append(filename)
  
  logger.info("")
  logger.info("Cleaning up and saving extra files to prep_files.tar")
  logger.debug("The files are: %s"%" ".join(tarlist2))
  subprocess.call("tar -cf prep_files.tar %s"%" ".join(tarlist2),shell=True)
  subprocess.call("rm -f %s"%" ".join(tarlist2),shell=True)

def _wizard(settings) :

  print "\nBecause you initiated protoms.py without any arguments, a wizard will walk you through the setup\n"
  
  print "What kind of simulation would you like to setup?"
  print "\t1) Equilibration"
  print "\t2) Sampling"
  print "\t3) Dual topology free energy"
  print "\t4) Single topology free energy"
  print "\t5) Grand-canonical Monte Carlo (GCMC)"
  print "\t6) JAWS stage 1 (Just Add Waters)"
  print "\t7) JAWS stage 2 (Just Add Waters)"
  print ">",
  valid = ["","1","2","3","4","5","6","7"]
  instr = raw_input()
  if instr in ["jon","dev"] :
    import matplotlib.pylab as plt
    img = np.load(simulationobjects.standard_filename(".ee.npz","tools"))[instr]
    plt.imshow(img)
    plt.show()
    return
  while instr not in valid :
    print "Please type a number between 1 and 7!"
    print ">",
    instr = raw_input()
  if instr == "" : return
  val = int(instr)
  vals = ["equilibration","sampling","dualtopology","singletopology","gcmc","jaws1","jaws2"]
  settings.simulation = vals[val-1]
  
  print "\nDo you have a protein that you would like to setup and simulate?"
  print "\tPlease enter the protein name or a PDB filename.\n\tPress enter if you don't have a protein."
  print ">",
  instr = raw_input()
  if instr != "" : settings.protein = instr

  print "\nDo you have a ligand that you would like to setup and simulate?"
  if settings.simulation not in ["dualtopology","singletopology"] :
    print "\tPlease enter the ligand name or a PDB filename. \n\tPress enter if you don't have a ligand."
    print ">",
    instr = raw_input()
    if instr != "" : args.ligand = [instr]
  else :
    while instr == "" :
      print "For the type of simulation you have selected, two ligands need to be given. \n\tPlease enter a ligand name or a PDB filename"
      print ">",
      instr = raw_input()
    args.ligand = [instr]     
    print "\tPlease enter the second ligand name or a PDB filename."
    instr = ""
    while instr == "" :
      print ">",
      instr = raw_input()
      if instr == "" :
        print "\tPlease note that absolute free energies cannot currently be set up via the wizard.\n\tPlease see protoms.py --fullhelp for details of how to set these up!"
        print "\tPlease enter the second ligand name or a PDB filename."
    args.ligand.append(instr)

  print "\nDo you have a co-factor or another solute that you would like to setup and simulate?"
  print "\tPlease enter the solute name or a PDB filename."
  print "\tPress enter if you don't have a solute or when you have specified all solutes."
  print ">",
  instr = raw_input()  
  if instr != "" :
    if args.ligand is None or len(args.ligand) == 0 :
      args.ligand = [instr]
    else :
      args.ligand.append(instr)
  while instr != "" : 
    print ">",
    instr = raw_input()  
    if instr != "" :
      args.ligand.append(instr)
    

if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = simulationobjects.MyArgumentParser(description="Program setup and run a ProtoMS simulations")
  parser.add_argument('-s','--simulation',choices=["none","equilibration","sampling","dualtopology","singletopology","gcmc","jaws1","jaws2"],help="the kind of simulation to setup",default="none")
  parser.add_argument('-f','--folders',nargs="+",help="folders to search for files ",default=["."])
  parser.add_argument('-p','--protein',help="the prefix of the protein")
  parser.add_argument('-l','--ligand',nargs="+",help="the prefix of the ligand(s)")
  parser.add_argument('-w','--water',help="the prefix of the water/solvent",default="water")
  parser.add_argument('-c','--cmdfile',help="the prefix of the command file",default="run")
  parser.add_argument('-sc','--scoop',help="the name of your protein scoop")
  parser.add_argument('-t','--template',nargs="+",help="the template files for your ligands")
  parser.add_argument('-r','--repeats',help="the number of repeats to be run (if more than 1) or a name for your repeat",default="") 
  # General control variables
  cntrlgroup = parser.add_argument_group("General control variables")
  cntrlgroup.add_argument('--outfolder',help="the ProtoMS output folder",default="")
  cntrlgroup.add_argument('--atomnames',help="a file with atom name conversions")
  cntrlgroup.add_argument('--watmodel',help="the name of the water model. Default = tip4p",choices=[ 'tip3p', 'tip4p'],default='tip4p')
  cntrlgroup.add_argument('--waterbox',help="a file with pre-equilibrated water molecules")
  cntrlgroup.add_argument('--setupseed',help="optional seed for random number generators in setup",default=None,type=int) 
  # Ligand setup variables
  liggroup = parser.add_argument_group("Ligand setup variables")
  liggroup.add_argument('--charge',nargs="+",type=float,help="the net charge of each ligand")
  liggroup.add_argument('--singlemap',help="the correspondance map for single-topology")
  liggroup.add_argument('--gaff',help="the version of GAFF to use for ligand", default="gaff16")
  # Protein setup variables
  protgroup = parser.add_argument_group("Protein setup variables")
  protgroup.add_argument('--center',help="the center of the scoop, if ligand is not available, either a string or a file with the coordinates",default=None)
  protgroup.add_argument('--innercut',type=float,help="maximum distance from ligand defining inner region of the scoop",default=16.0)
  protgroup.add_argument('--outercut',type=float,help="maximum distance from ligand defining outer region of the scoop",default=20.0)
  protgroup.add_argument('--flexin',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the inner region",default="flexible")
  protgroup.add_argument('--flexout',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the outer region",default="sidechain")
  protgroup.add_argument('--scooplimit',help="the minimum difference between number of residues in protein and scoop for scoop to be retained",default=10)
  protgroup.add_argument('--capradius',type=float,help="the radius of the droplet around the protein",default=30.0)
  # Simulation parameters
  simgroup = parser.add_argument_group("Simulatiom parameters")
  simgroup.add_argument('--lambdas',nargs="+",type=float,help="the lambda values or the number of lambdas",default=[16])
  simgroup.add_argument('--adams',nargs="+",type=float,help="the Adam/B values for the GCMC",default=0)
  simgroup.add_argument('--adamsrange',nargs="+",type=float,help="the upper and lower Adam/B values for the GCMC and, optionally, the number of values desired (default value every 1.0), e.g. -1 -16 gives all integers between and including -1 and -16",default=None)
  simgroup.add_argument('--gcmcwater',help="a pdb file with a box of water to do GCMC on or an integer corresponding to the number of water molecules to add")
  simgroup.add_argument('--gcmcbox',nargs="+",help="a pdb file with box dimensions for the GCMC box, or a list of origin(x,y,z) and length(x,y,z) coordinates")
  simgroup.add_argument('--jawsbias',type=float,nargs="+",help="the bias in JAWS-2",default=[6.5])
  simgroup.add_argument('--nequil',type=float,help="the number of equilibration steps",default=5E6)
  simgroup.add_argument('--nprod',type=float,help="the number of production steps",default=40E6)
  simgroup.add_argument('--dumpfreq',type=float,help="the output dump frequency",default=1E5)
  simgroup.add_argument('--ranseed',help="the value of the random seed you wish to simulate with. If None, then a seed is randomly generated. Default=None",default=None)
  simgroup.add_argument('--absolute',action='store_true',help="whether an absolute free energy calculation is to be run. Default=False",default=False)
  simgroup.add_argument('--dovacuum',action='store_true',help="turn on vacuum simulation for simulation types equilibration and sampling",default=False)
  simgroup.add_argument('--testrun',action='store_true',help="setup a short test run. Default=False",default=False)
  simgroup.add_argument('--cleanup',action='store_true',help="Clean up extra files. Default=False",default=False)
  simgroup.add_argument('--tune',action='store_true',help='Carry out dihedral tuning simulation',default=False)
  simgroup.add_argument('--softcore', type=str, default='auto', 
                        choices=('auto', 'all', 'none', 'manual'),
                        help="determine which atoms to apply softcore potentials to.\n "
                             "'all'=softcores applied to all atoms of both solutes, "
                             "'none'=softcores not applied to any atoms\n "
                             "'mixed'=softcores will be applied only to non matching "
                             "atoms within ligand structures")
  parser.add_argument(
    '--spec-softcore', type=str,
    help='Specify atoms to add or remove from softcore selections. Can be '
         'up to two, space separated, strings of the form "N:AT1,AT2,-AT3". '
         'N should be either "1" or "2" indicating the corresponding ligand. '
         'The comma separated list of atom names are added to the softcore '
         'selection. A preceding dash for an atom name specifies it should be'
         ' removed from the softcore selection.')
  args = parser.parse_args()
 
  print r"""
            __|\
         .-'    '-.
        / .--, _ a L
      .J (  '-' "'--'
     '-'-.)  .~~~~~~~~~~~~~~~~~~~~.
             |                    |     __
             |    Welcome to      | ,.-'e ''-'7
             |    protoms.py      |  '--.    (
             |                    |      ),   \
             '~~~~~~~~~~~~~~~~~~~~'      ` )  :
                                      ,__.'_.'
                                      '-, (
                                        '--'  """
  print ""

  # Setup the logger - logger is global
  logger.debug("Running protoms.py at %s"%time.strftime("%d/%m/%Y - %H:%M:%S"))
  logger.debug("Command line arguments = %s"%" ".join(sys.argv[1:]))
  logger.debug("Settings = %s"%args)

  np.random.seed(args.setupseed)

  # Adds current folder to the folders
  args.folders.append(".")
  
  # If we run without any command-line arguments, initiate the wizard
  if len(sys.argv) == 1 :
    _wizard(args)

  if args.protein is None and args.ligand is None and args.scoop is None:
    print "Nothing to do, so exit!"
    quit()
 
  # Set $PROTOMSHOME 
  if os.getenv("PROTOMSHOME") is None :
    string = os.path.dirname(os.path.abspath(__file__))
    logger.info("Setting PROTOMSHOME to %s"%string)
    os.environ["PROTOMSHOME"] = string # This does not change the original shell

  # Try to find a default water box
  if args.waterbox is None :
    args.waterbox = simulationobjects.standard_filename("wbox_"+args.watmodel.lower()+".pdb","data")
  if not os.path.isfile(args.waterbox) : 
    msg = "Could not find file (%s) with pre-equilibrated waters"%args.waterbox
    logger.error(msg)
    raise simulationobjects.SetupError(msg)

  # Setup list of files to be stored away
  tarlist = []
     
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

    # Make a dummy PDB structure
    if args.simulation == "dualtopology" and len(ligands) < 2 :
      if args.absolute :
#    if args.simulation == "dualtopology" and args.absolute :
        ligands.insert(1,"*dummy")
        prefix0 = ligands[0]
        dummy_name = os.path.basename(_get_prefix(prefix0))+"_dummy.pdb"
        ligand_files["*dummy"] = {}
        ligand_files["*dummy"]["pdb"] = dummy_name
        ligand_files["*dummy"]["obj"] = tools.make_dummy(ligand_files[prefix0]["obj"])
        ligand_files["*dummy"]["obj"].write(ligand_files["*dummy"]["pdb"])
        ligand_files["*dummy"]["tem"] = simulationobjects.standard_filename("dummy.tem","data")
        logger.info("")
        logger.info("Creating dummy PDB-file for ligand: %s"%ligand_files["*dummy"]["pdb"])
      else :
        msg = "You are trying to run a dual topology simulation but with only one ligand! (%s.pdb)\nIf you meant to run an absolute free energy calculation, rerun with the --absolute option."%ligands[0]
        logger.error(msg)
        raise simulationobjects.SetupError(msg)


    # Create merged pdb objects
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
      if l[0] == "*" : continue # Skip ligands created in the script, i.e. the dummy
      if args.charge is not None and i < len(args.charge): 
        charge = args.charge[i]
      else :
        charge = 0
      if i > 1 : ligobj12 = None
      prefix = _get_prefix(l)
      _prep_ligand(ligand_files[prefix],i==0,charge,ligobj12,args.folders,tarlist,args) 

    ligpdbs = [ligand_files[l]["pdb"] for l in ligands]
    ligtems = [ligand_files[l]["tem"] for l in ligands]
    ligand_water = ligand_files[ligands[0]]["wat"]
 
    # Here we need to make single topology templates, if requested
    if args.simulation == "singletopology" :
      ligtems,ligtems2,ligtems3 = _prep_singletopology(ligpdbs,ligtems,tarlist,args)

    # Here we will merge ligand template files if there is more than one
    if len(ligtems) > 1 :
      logger.info("")
      ligtems = _merge_templates(ligtems,tarlist)
      if args.simulation == "singletopology" : 
        ligtems2 = _merge_templates(ligtems2,tarlist)    
        ligtems3 = _merge_templates(ligtems3,tarlist)    
    
  # Prepare the protein
  protein_file = None
  water_file = None
  if args.protein is not None or args.scoop is not None:
    protein_file,water_file = _prep_protein(args.protein,ligobjs,args.water,args.folders,tarlist,args)

  # Extra preparation for GCMC or JAWS-1     
  if args.simulation in ["gcmc","jaws1"] :
    if water_file is None:
        msg = "GCMC and JAWS1 not supported without protein or scoop"
        logger.error(msg)
        raise simulationobjects.SetupError(msg)
    args.gcmcwater,water_file = _prep_gcmc(ligands,ligand_files,water_file,tarlist,args)
  
  # Extra preparation for JAWS-2
  if args.simulation == "jaws2"  :
    single_wat,other_wat,water_file = _prep_jaws2(water_file,tarlist,args)

  # Check of test run
  if args.testrun :
    if args.nequil == 5E6 :
      args.nequil = 0
    if args.nprod == 40E6 :
      args.nprod = 4000
    if args.dumpfreq == 1E5 :
      args.dumpfreq = 10
    if len(args.lambdas) == 1 and args.lambdas[0] == 16 :
      args.lambdas = [4]
      
  # Create ProtoMS command files
  ranseed=args.ranseed
  if args.simulation == "singletopology" :
   postfix = ["_ele","_vdw","_comb"]
  elif args.simulation == "jaws2" :
    postfix = ["_jaws2-w%d"%(i+1) for i in range(len(single_wat.pdbs))] 
    if args.outfolder == "" : args.outfolder = "out"
  else :
   postfix = [""]
  if args.repeats.isdigit():
    args.repeats = range(1,int(args.repeats)+1)
  else :
    args.repeats = [args.repeats.lower()]

  repeats = []
  for post in postfix :
    for repeat in args.repeats:
      repeats.append(str(repeat) + post)

  outfolder = args.outfolder
  if args.outfolder == "" : 
    if args.simulation == "gcmc" :
      outfolder = "out_gcmc"
    elif args.simulation in ["jaws1","jaws2"] :
      outfolder = "out_jaws"
    else :
      outfolder = "out"

  for repeat in repeats :
    args.outfolder = outfolder + repeat
    #setattr(args,"outfolder","out"+repeat)
    if not args.simulation in ["singletopology","jaws2"] or "_ele" in repeat :
      free_cmd,bnd_cmd,gas_cmd = tools.generate_input(protein_file,ligpdbs,ligtems,water_file,ligand_water,ranseed,args)
    elif args.simulation == "singletopology" and "_vdw" in repeat :
      free_cmd,bnd_cmd,gas_cmd = tools.generate_input(protein_file,ligpdbs,ligtems2,water_file,ligand_water,ranseed,args)
    elif args.simulation == "singletopology" and "_comb" in repeat :
      free_cmd,bnd_cmd,gas_cmd = tools.generate_input(protein_file,ligpdbs,ligtems3,water_file,ligand_water,ranseed,args)
    elif args.simulation == "jaws2" :
      idx = int(repeat.split("-")[-1][1:])
      args.gcmcwater = "jaws2_wat%d.pdb"%idx
      jaws2wat = "jaws2_not%d.pdb"%idx
      free_cmd,bnd_cmd,gas_cmd = tools.generate_input(protein_file,ligpdbs,ligtems,water_file+" "+jaws2wat,ligand_water,ranseed,args)

    args.cmdfile=args.cmdfile.lower() #protoMS cannot handle cmd files containing upper case letters
    if free_cmd is not None : 
      free_cmd.writeCommandFile(args.cmdfile+repeat+"_free.cmd")
    if bnd_cmd is not None :
      if args.simulation == "gcmc" :
        bnd_cmd.writeCommandFile(args.cmdfile+repeat+"_gcmc.cmd") 
      elif args.simulation in ["jaws1","jaws2"] :
        bnd_cmd.writeCommandFile(args.cmdfile+repeat+"_jaws.cmd") 
      else :
        bnd_cmd.writeCommandFile(args.cmdfile+repeat+"_bnd.cmd")       
    if gas_cmd is not None :
      if args.absolute:
        # in this case, gas_cmd contains a cmd file to account for introduction of the harmonic restraint
        gas_cmd.writeCommandFile(args.cmdfile+repeat+"_bnd_rstr.cmd")   
      else:
        gas_cmd.writeCommandFile(args.cmdfile+repeat+"_gas.cmd")   
      
    
  if args.cleanup :
    _cleanup(tarlist)
    


