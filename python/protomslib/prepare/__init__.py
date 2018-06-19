import distutils.spawn
import logging
import numpy as np
import os
import six
import subprocess
import tempfile
from .. import simulationobjects
from ..gcmc import clear_gcmcbox, distribute_particles
from ..gcmc import make_gcmcbox, print_bequil
from ..solvate import solvate
from ..templates import merge_templates, build_template
from ..templates import make_single, summarize_single, write_map
from ..utils import _get_prefix, _locate_file
from .scoop import scoop
from .water import convertwater, split_waters, set_jaws2_box

logger = logging.getLogger('protoms')


def _merge_templates(templates, tarlist):
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
    for tem in templates:
        tarlist.append(
            tem)  # All of the original templates can safely be stored away
    temfile = merge_templates(templates)
    if isinstance(temfile, simulationobjects.TemplateFile):
        allnames = "-".join(t.name.lower() for t in temfile.templates)
        templates = [allnames + ".tem"]
        temfile.write(allnames + ".tem")
        logger.info("Created a re-numbered template files for all ligands: %s"
                    % templates[0])
        return templates
    else:
        return [temfile]


def _load_ligand_pdb(ligprefix, folders):
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
    pdbfile = _locate_file(ligprefix + ".pdb", folders)

    # Cannot do anything without a pdb-file so raise an exception
    if pdbfile is None:
        msg = "Ligand file %s.pdb could not be found" % ligprefix
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    return pdbfile, simulationobjects.PDBFile(filename=pdbfile)


def _prep_ligand(files, first, charge, ligobj12, folders, tarlist, settings):
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
    files["prepi"] = _locate_file(ligprefix + ".prepi", folders)
    files["frcmod"] = _locate_file(ligprefix + ".frcmod", folders)
    files["zmat"] = _locate_file(ligprefix + ".zmat", folders)
    files["wat"] = _locate_file(ligprefix + "_box.pdb", folders)

    logger.info("")
    logger.info("Setting up ligand: %s..." % files["pdb"])

    if len(files["obj"].residues) < 1:
        if len(files["obj"].solvents) > 0:
            raise simulationobjects.SetupError(
                "Your ligand in %s is recognized as solvent. "
                "Please change the residue name." % files["obj"].name)
        else:
            raise simulationobjects.SetupError(
                "No residues found in %s." % files["obj"].name)

    # Get the ligand name from the pdb header
    if "HEADER" in files["obj"].header:
        words = files["obj"].header.strip().split()
        ligname = words[words.index("HEADER") + 1]
    else:
        logger.info(
            "Unable to find header information in the PDB-file, "
            "adding it automatically."
        )
        ligname = files["obj"].residues[1].name
        files["obj"].header = files["obj"].header + "HEADER " + ligname + "\n"
        files["obj"].write(files["pdb"])

    # Try to locate template for the ligand
    files["tem"] = None
    if settings.template is None:
        files["tem"] = _locate_file(ligprefix + ".tem", folders)
    else:
        for f in settings.template:
            tempfile = _locate_file(f, folders)
            if tempfile is None:
                files["tem"] = tempfile
            else:
                tem = simulationobjects.TemplateFile(filename=tempfile)
                if ligname.lower() in [t.name.lower() for t in tem.templates]:
                    files["tem"] = tempfile
                    break

    # Check to see if we have a template file
    if files["tem"] is None:
        resnam = files["obj"].residues[
            1].name  # Set the prepi name and template name to the residue name
        if files["prepi"] is None:
            # Here we need to run Antechamber
            logger.info(
                "Running antechamber. Please check the output carefully")
            files["prepi"] = run_antechamber(files["pdb"], charge, resnam)
            logger.info("Created prepi-file: %s" % files["prepi"])
            tarlist.append(files["prepi"])
        if files["frcmod"] is None:
            # Here we need to run parmchk
            logger.info("Running parmchk. Please check the output carefully")
            files["frcmod"] = run_parmchk(files["pdb"])
            logger.info("Created frcmod-file: %s" % files["frcmod"])
            tarlist.append(files["frcmod"])

        # By this stage we have all necessary files to make the template file
        files["tem"] = ligprefix + ".tem"
        tem = build_template(
            temfile=files["tem"],
            prepifile=files["prepi"],
            zmatfile=files["zmat"],
            frcmodfile=files["frcmod"],
            resname=resnam,
            gaffversion=settings.gaff)
        tem.write(files["tem"])
        if files["zmat"] is None:
            files["zmat"] = ligprefix + ".zmat"
            tem.templates[0].write_zmat(files["zmat"])
            logger.info(
                "Created zmatrix (%s) for ligand. Please check the output"
                " carefully" % files["zmat"])
            tarlist.append(files["zmat"])
        logger.info(
            "Created ProtoMS template-file (%s) for ligand."
            " Please check the output carefully" % files["tem"])

    # Check to see if we have solvated the ligand
    if files["wat"] is None:

        # Setting the solute
        if ligobj12 is None:
            solute = files["obj"]
        else:
            solute = ligobj12
        # Calling the routine
        files["wat"] = ligprefix + "_box.pdb"
        boxpdb = solvate(
            settings.waterbox,
            ligand=solute,
            protein=None,
            geometry="box",
            padding=10.0,
            radius=30.0,
            center="cent",
            namescheme="ProtoMS")

        boxpdb.write(files["wat"])
        logger.info("Created waterbox-file: %s" % (files["wat"]))
        if (settings.simulation in ["dualtopology", "singletopology"]
                and not first) or (settings.simulation in [
                    "gcmc", "jaws1", "jaws2"
                ]):
            tarlist.append(files["wat"])

    return files


def _prep_protein(protprefix, ligands, watprefix, folders, tarlist, settings):
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

    if protprefix is not None:
        protprefix = _get_prefix(protprefix)
    else:
        protprefix = _get_prefix(settings.scoop)
    watprefix = _get_prefix(watprefix)
    if settings.scoop is None:
        scoopprefix = protprefix + "_scoop"
    else:
        scoopprefix = _get_prefix(settings.scoop)

    # Try to locate necessary protein files
    protein_orig_file = _locate_file(protprefix + ".pdb", folders)
    protein_pms_file = _locate_file(protprefix + "_pms.pdb", folders)
    protein_scoop_file = _locate_file(scoopprefix + ".pdb", folders)
    protein_water = _locate_file(watprefix + ".pdb", folders)

    # Cannot do anything without a pdb-file, so raise an exception
    if settings.protein is None and protein_scoop_file is None:
        msg = "Specified scoop (%s.pdb) could not be found" % scoopprefix
        logger.error(msg)
        raise simulationobjects.SetupError(msg)
    if protein_orig_file is None and protein_scoop_file is None and \
       protein_pms_file is None:
        msg = "Protein file (%s.pdb) and protein scoop file (%s.pdb) and " \
              "protein pms file (%s_pms.pdb) could not be found" % (
                  protprefix, scoopprefix, protprefix)
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    logger.info("")
    logger.info("Setting up protein: %s..." % protein_orig_file)

    if settings.scoop is not None and protein_scoop_file is None:
        logger.info(
            "Specified scoop (%s.pdb) not found. Ignoring..." % scoopprefix)

    protobj = None
    if protein_scoop_file is None and protein_pms_file is None:
        # Start with an object for the original pdb-file
        protobj = simulationobjects.PDBFile(filename=protein_orig_file)

        # Trying to locate a default conversions file
        if settings.atomnames is None:
            conversionfile = simulationobjects.standard_filename(
                "atomnamesmap.dat", "data")
        else:
            conversionfile = settings.atomnames
        if not os.path.isfile(conversionfile):
            msg = "Could not find file (%s) with atom name conversions" % \
                conversionfile
            logger.error(msg)
            raise simulationobjects.SetupError(msg)

        # Converting to ProtoMS atom names
        protobj2 = pdb2pms(protobj, "amber", conversionfile)
        abortedconv = protobj2 == protobj  # Indicates abortion of conversion
        protobj = protobj2

        # Converting water molecules to specified model
        protobj = convertwater(protobj, settings.watmodel)

        # Defining the center of the scoop...
        if ligands is None and settings.gcmcbox is None:
            if settings.center is not None:
                ligobj = settings.center
            else:
                protobj.getCenter()
                ligobj = "%f %f %f" % tuple(protobj.center)
                logger.warning(
                    "Warning: No specified center for protein scoop. "
                    "Using the center of the protein."
                )
        elif ligands is None and settings.gcmcbox is not None:
            ligobj = simulationobjects.PDBFile(filename=settings.gcmcbox[0])
            logger.info("Scooping protein around %s..." % settings.gcmcbox[0])
        else:
            ligobj = ligands

        # Here we need to call the routine to make a scoop
        protobj_scooped = scoop(
            protobj,
            ligobj,
            innercut=settings.innercut,
            outercut=settings.outercut,
            flexin=settings.flexin,
            flexout=settings.flexout,
            scooplimit=settings.scooplimit)

        nresdiff = len(protobj.residues) - len(protobj_scooped.residues)
        protein_scoop_file = _get_prefix(protein_orig_file) + "_scoop.pdb"
        protobj = protobj_scooped
        if nresdiff < settings.scooplimit:
            if not abortedconv:
                protein_pms_file = _get_prefix(protein_orig_file) + "_pms.pdb"
                logger.info("Created %s instead." % protein_pms_file)
                header = protobj.header
                header += "REMARK Atoms renamed according to ProtoMS naming" \
                    " standards.\n"
                protobj.write(
                    filename=protein_pms_file, header=header, solvents=False)
            else:
                protein_pms_file = protein_orig_file
            tarlist.append(protein_scoop_file)
            protein_scoop_file = None
        else:
            logger.info("Created scoop-pdb file by removing %d residues: %s" %
                        (nresdiff, protein_scoop_file))
            protobj.write(protein_scoop_file, renumber=True, solvents=False)

    if protein_water is None:

        if protobj is None:
            if protein_scoop_file is not None:
                solute = protein_scoop_file
            else:
                solute = protein_pms_file
        else:
            solute = protobj

        # Calling the routine
        protein_water = watprefix + ".pdb"
        cappdb = solvate(
            settings.waterbox,
            ligand=ligands,
            protein=solute,
            geometry="droplet",
            padding=10.0,
            radius=settings.capradius,
            center="cent",
            namescheme="ProtoMS")
        cappdb.write(protein_water)
        logger.info("Created water cap-file: %s" % (protein_water))

    if protein_pms_file is not None:
        return protein_pms_file, protein_water
    else:
        return protein_scoop_file, protein_water


def _prep_singletopology(pdbs, templates1, tarlist, settings):
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
        msg = "Single topology calculations require two ligand templates to " \
              "perturb between.\nYou only have one: %s \nIf you are trying " \
              "to run an absolute free energy calculation\n(i.e. perturb to" \
              " nothing), please use dual topology." % templates1[0]
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    # The original templates file can now be stored away
    tarlist.append(templates1[0])
    tarlist.append(templates1[1])

    logger.info("")
    logger.info(
        "Setting up single-topology correspondance map and templates...")
    eletem, vdwtem, combtem, cmap = make_single(
        tem1,
        tem2,
        pdbs[0],
        pdbs[1],
        settings.singlemap,
        gaffversion=settings.gaff)
    summarize_single(eletem, vdwtem, logger.debug)

    pdbs.pop(1)
    templates1.pop(1)
    templates2 = list(templates1)
    templates3 = list(templates1)

    prefix = os.path.join(
        os.path.dirname(pdbs[0]),
        tem1.templates[0].name.lower() + "t" + tem2.templates[0].name.lower())
    templates1[0] = prefix + "_ele.tem"
    templates2[0] = prefix + "_vdw.tem"
    templates3[0] = prefix + "_comb.tem"

    eletem.write(templates1[0])
    vdwtem.write(templates2[0])
    combtem.write(templates3[0])
    logger.info("")
    logger.info("Created template %s for electrostatic perturbation. "
                "Please check the output carefully." % templates1[0])
    logger.info("Created template %s for van der Waals-perturbation. "
                "Please check the output carefully." % templates2[0])
    logger.info("Created template %s for combined perturbation. "
                "Please check the output carefully." % templates3[0])

    if settings.singlemap is None:
        settings.singlemap = "single_cmap.dat"
        tarlist.append(settings.singlemap)
    write_map(cmap, settings.singlemap)
    logger.info("Saved correspondance map to: %s" % settings.singlemap)

    return templates1, templates2, templates3


def _prep_gcmc(ligands, ligand_files, waters, tarlist, settings):
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

    def pdb2box(pdbobj, padding=2.0):
        boxpdb = "%s_box.pdb" % settings.simulation
        box = make_gcmcbox(pdbobj, boxpdb, padding)
        simulationobjects.write_box(boxpdb, box)
        logger.info("")
        logger.info("Created %s to visualize GCMC/JAWS-1 simulation box. "
                    "Please check the output carefully" % boxpdb)
        return boxpdb

    def arrange_wats(box, water_file, out_name):
        # find smallest box including all oxygens in my water_file
        waterobj = simulationobjects.PDBFile(filename=water_file)
        if len(waterobj.residues) != 0:
            waterbox = waterobj.getBox(
                atomlist=[next(six.itervalues(waterobj.residues)).atoms[0].name])
        else:
            waterbox = waterobj.getBox(
                atomlist=[next(six.itervalues(waterobj.solvents)).atoms[0].name])

        # check whether the box of the water oxygens is
        # within the box provided as argument
        box_origen_below = np.all(box['origin'] < waterbox['origin'])
        box_end_avobe = np.all(
            box['origin'] + box['len'] > waterbox['origin'] + waterbox['len'])
        if not (box_origen_below and box_end_avobe) or (
                np.all(waterbox['len'] < 0.5) and len(waterobj.solvents) > 1):
            # re-distribute the waters if required according to provided box
            logger.info("\nRedistributing your GCMC / JAWS1 waters"
                        " according to the specified box\n")
            arranged_obj = distribute_particles(box, water_file)
            arranged_obj.write(filename=out_name)
            return arranged_obj, out_name
        else:
            # make sure the header in the waterobj is the box dimensions
            # which will be used in this simulation
            box_extremes = tuple(box['origin']) + tuple(
                np.add(box['origin'], box['len']))
            waterobj.header = "HEADER box %.3f %.3f %.3f %.3f %.3f %.3f \n" \
                % box_extremes
            return waterobj, out_name

    def fill_box(settings, boxpdb, box, ghost_name):
        # Fill the box with waters
        if settings.gcmcwater is None:
            # Creating ghost waters with a density 1.5 times that of bulk water
            box_volume = box['len'][0] * box['len'][1] * box['len'][2]
            ghostnum = str(
                int(np.ceil(box_volume * 0.0334 * 1.5))
            )  # 0.0334 waters per Angs.^3 is the number density of bulk water.
            ghostobj = distribute_particles(
                box=box, particles=ghostnum, watermodel=settings.watmodel)
            ghostobj.write(filename=ghost_name)
        elif settings.gcmcwater.isdigit():
            ghostobj = distribute_particles(
                box, settings.gcmcwater, watermodel=settings.watmodel)
            ghostobj.write(filename=ghost_name)
        else:
            ghostobj, ghost_name = arrange_wats(box, settings.gcmcwater,
                                                ghost_name)
        return ghostobj, ghost_name

    # Check consistency of command-line arguments
    if settings.gcmcwater is not None and not settings.gcmcwater.isdigit():
        gcmcwater = _locate_file(settings.gcmcwater, settings.folders)
        if gcmcwater is None:
            msg = "File %s given as gcmcwater could not be found." \
                % settings.gcmcwater
            logger.error(msg)
            raise simulationobjects.SetupError(msg)
        settings.gcmcwater = gcmcwater
    if settings.gcmcbox is not None and len(settings.gcmcbox) is 1:
        gcmcbox = _locate_file(settings.gcmcbox[0], settings.folders)
        if gcmcbox is None:
            msg = "File %s given as gcmcbox could not be found." \
                % settings.gcmcbox[0]
            logger.error(msg)
            raise simulationobjects.SetupError(msg)
        settings.gcmcbox = gcmcbox
    elif settings.gcmcbox is not None and len(settings.gcmcbox) < 6:
        msg = "6 arguments expected to define the GCMC/JAWS1 box dimensions" \
            ", %d provided: %s" % (
                len(settings.gcmcbox), " ".join(settings.gcmcbox))
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    ghost_name = "%s_wat.pdb" % settings.simulation
    write = True

    # If the gcmcbox has been set as a file
    if settings.gcmcbox is not None and isinstance(settings.gcmcbox, str):
        gcmcboxobj = simulationobjects.PDBFile(filename=settings.gcmcbox)
        # Try to find the box in the header of the file
        box = simulationobjects.find_box(gcmcboxobj)
        if box is None:
            boxpdb = pdb2box(gcmcboxobj, padding=0.0)
            gcmcboxobj = simulationobjects.PDBFile(filename=boxpdb)
            box = simulationobjects.find_box(gcmcboxobj)
        else:
            boxpdb = settings.gcmcbox
        if "center" in box:
            box['origin'] = np.array([
                coord - box["len"][ind] / 2
                for ind, coord in enumerate(box["center"])
            ])

        # Fill the box with waters
        ghostobj, ghost_name = fill_box(settings, boxpdb, box, ghost_name)

    # If the dimensions of the gcmc box have been provided
    elif settings.gcmcbox is not None and len(settings.gcmcbox) is 6 and [
            _is_float(value) for value in settings.gcmcbox
    ]:
        # Generate the box pdb
        boxpdb = "%s_box.pdb" % settings.simulation
        try:
            box = {}
            box['origin'] = np.array(
                [float(value) for value in settings.gcmcbox[:3]])
            box['len'] = np.array(
                [float(value) for value in settings.gcmcbox[3:]])
        except Exception:
            raise simulationobjects.SetupError(
                "The box dimensions %s could not be correctly interpreted" %
                settings.gcmcbox)
        simulationobjects.write_box(boxpdb, box)
        logger.info("")
        logger.info("Created %s to visualize GCMC/JAWS-1 simulation box. "
                    "Please check the output carefully" % boxpdb)
        # Fill the box with waters
        ghostobj, ghost_name = fill_box(settings, boxpdb, box, ghost_name)

    # If no box information has been provided but we have a gcmcwater file
    elif (settings.gcmcbox is None and settings.gcmcwater is not None
          and not settings.gcmcwater.isdigit()):
        waterobj = simulationobjects.PDBFile(filename=settings.gcmcwater)
        box = simulationobjects.find_box(waterobj)
        if box is None:
            raise simulationobjects.SetupError(
                "Cannot set up simulation without a box.")
        else:
            ghostobj = waterobj
            ghost_name = waterobj.name
            write = False

    # If no information on the box has been provided
    elif ligands:
        # Use the ligand to define the GCMC box
        boxpdb = pdb2box(ligand_files[ligands[0]]["obj"])
        gcmcboxobj = simulationobjects.PDBFile(filename=boxpdb)
        box = simulationobjects.find_box(gcmcboxobj)
        if "center" in box:
            box['origin'] = np.array([
                coord - box["len"][ind] / 2
                for ind, coord in enumerate(box["center"])
            ])
        # Fill the box with waters
        ghostobj, ghost_name = fill_box(settings, boxpdb, box, ghost_name)
    else:
        msg = "Cannot define a GCMC or JAWS-1 simulation box without a" \
              " ligand and without the gcmcbox setting"
        logger.error(msg)
        simulationobjects.SetupError(msg)

    print_bequil(box["len"])
    # Write the waters to disc
    if write:
        ghostobj.write(ghost_name)
        logger.info("")
        logger.info("Created %s; it contains the GCMC or JAWS-1 simulation "
                    "waters. Please check the output carefully" % ghost_name)

    # Clear the GCMC/JAWS-1 box from solvation waters
    logger.info("")
    nrem, waters2 = clear_gcmcbox(ghostobj, waters)
    if nrem > 0:
        waters2_name = _get_prefix(str(waters2)) + "_clr.pdb"
        logger.info("Created water cap-file: %s" % waters2_name)
        tarlist.append(waters)
        waters2.write(waters2_name)
        waters = waters2_name

    return ghost_name, waters


def _prep_jaws2(water_file, tarlist, settings):
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

    if settings.gcmcwater is None:
        msg = "You must set gcmcwater settings when preparing JAWS-2 input"
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    single_wat, other_wat = split_waters(settings.gcmcwater)
    single_wat = set_jaws2_box(single_wat)
    logger.info("")
    logger.info("Creating water PDB-files for JAWS-2 called "
                "jaws2_wat*.pdb and jaws2_not*.pdb")
    single_wat.write(
        ["jaws2_wat%d.pdb" % (i + 1) for i in range(len(single_wat.pdbs))])
    other_wat.write(
        ["jaws2_not%d.pdb" % (i + 1) for i in range(len(single_wat.pdbs))])

    nrem = 0
    for count, watobj in enumerate(single_wat.pdbs):
        # Clear the JAWS-2 box from solvation waters
        watobj.name = "jaws2_wat%d.pdb" % (count + 1)
        n, water_file = clear_gcmcbox(watobj, water_file)
        nrem = nrem + n
    if nrem > 0:
        waters2_name = _get_prefix(str(water_file)) + "_clr.pdb"
        logger.info("Created water cap-file: %s" % waters2_name)
        tarlist.append(water_file)
        water_file.write(waters2_name)
        water_file = waters2_name
    else:
        water_file = water_file.name

    return single_wat, other_wat, water_file


def _run_program(name, command):
    """
    Wrapper for an AmberTools program with default parameters

    run_program ( name, command )

    Parameters
    ----------
    name : string
      the name of the program to run
    command : string
      the command to execute it

    Raises
    ------
    SetupError
      if the program failed to execute properly
    """

    # Create temporary file to write stdout, stderr of the program
    tmpfile, tmpname = tempfile.mkstemp()
    ret_code = subprocess.call(
        command, shell=True, stdout=tmpfile, stderr=tmpfile)
    # Catch some error codes
    if ret_code == 127:
        msg = "Unable to find executable, please make sure this is" \
            " present in your PATH or $AMBERHOME/bin."
        logger.error(msg)
        raise simulationobjects.SetupError(msg)
    if ret_code == 1:
        # Get the error message from the temporary file
        errmsg = "\n".join(line for line in open(tmpname).readlines())
        os.remove(tmpname)
        msg = "%s was not able to run successfully. Please check output. " \
            " Error message was:\n%s" % (name, errmsg)
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    os.remove(tmpname)


def _get_executable_path(name):
    """Robust routine that looks for Amber Tools executables first in
    $AMBERHOME/bin and then in the path generally.

    Parameters
    ----------
    name : string
      the name of the program to find

    Returns
    -------
    string
        the full path of the desired executable

    Raises
    ------
    SetupError
        if the executable could not be found
    """
    try:
        os.environ['AMBERHOME']
    except KeyError:
        raise simulationobjects.SetupError(
            "The environmental variable $AMBERHOME is not set. This is "
            "required to use components of AMBER Tools."
        )

    try:
        exe_path = os.environ['AMBERHOME'] + '/bin/' + name
        if not os.path.isfile(exe_path):
            raise KeyError
        return exe_path
    except KeyError:
        logger.debug(
            "Unable to find executable in $AMBERHOME/bin. Looking in PATH.")

    exe_path = distutils.spawn.find_executable(name)
    if exe_path is None:
        msg = "Unable to find %s executable, please make sure this is" \
            " present in your PATH or $AMBERHOME/bin." % name
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    return exe_path


def run_antechamber(lig, charge, resnam=None):
    """
    Wrapper for antechamber with default parameters and AM1-BCC charges

    Parameters
    ----------
    lig : PDBFile or string
      the ligand to run Antechamber on
    charge : int
      the net charge of the ligand
    resnam : string, optional
      the residue name of the ligand

    Returns
    -------
    string
      the filename of the created prepi file
    """

    logger.debug("Running run_antechamber with arguments: ")
    logger.debug("\tlig    = %s" % lig)
    logger.debug("\tcharge = %d" % charge)
    logger.debug("\tresnam = %s" % resnam)
    logger.debug(
        "This will generate an Amber prepi file with AM1-BCC and GAFF atom types"
    )

    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig, six.string_types):
        name = lig
    else:
        name = lig.name

    if resnam is None:
        resnamstr = ""
    else:
        resnamstr = "-rn " + resnam

    ante_exe = _get_executable_path('antechamber')

    # Remove the extension from the filename
    out_name = os.path.splitext(name)[0]
    cmd = '%s -i %s -fi pdb -o %s.prepi -fo prepi -c bcc -nc %d %s -pf y' % (
        ante_exe, name, out_name, charge, resnamstr)
    _run_program('antechamber', cmd)
    subprocess.call("rm sqm.in sqm.out sqm.pdb", shell=True)
    return "%s.prepi" % out_name


def run_parmchk(lig):
    """
    Wrapper for parmcheck with default parameters

    Parameters
    ----------
    lig : PDBFile or string
      the ligand to run Parmcheck on

    Returns
    -------
    string
      the filename of the created frcmod file
    """

    logger.debug("Running run_parmchk with arguments: ")
    logger.debug("\tlig = %s" % lig)
    logger.debug(
        "This will generate an Amber frcmod file with additional parameters")

    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig, six.string_types):
        name = lig
    else:
        name = lig.name

    parm_exe = _get_executable_path('parmchk')

    # Remove the extension from the filename
    out_name = os.path.splitext(name)[0]
    cmd = '%s -i %s.prepi -f prepi -o %s.frcmod' % (parm_exe, out_name,
                                                    out_name)
    _run_program('parmchk', cmd)
    return "%s.frcmod" % out_name


def read_convfile(file=None, inmode=None, outmode=None):
    """
    Reads a conversion file into a dictionary

    Parameters
    ----------
    file : string
      the filename of the conversion file
    inmode : string
      the style of the input PDB-file
    outmode : string
      the style of the output PDB-file

    Returns
    -------
    a dictionary of the conversion

    Raises
    ------
    ValueError
      if any of the parameters are none
    """

    if file is None:
        raise ValueError("You must pass file!")
    if inmode is None:
        raise ValueError("You must specify a mode!")
    if outmode is None:
        raise ValueError("You must specify a mode!")
    conv = {}
    stream = open(file, 'r')
    buffer = stream.readlines()
    stream.close()
    for line in buffer:
        if line.startswith('backbone'):
            key = 'backbone'
            conv[key] = {}
            continue
        if line.startswith('residue'):
            elems = line.split()
            key = elems[1]
            conv[key] = {}
            continue
        if line.startswith('atom'):
            elems = line.split()
            # Locate inmode
            found = False
            for x in range(0, len(elems)):
                if elems[x] == inmode:
                    found = True
                    old = elems[x + 1]
                    break
            if not found:
                continue
            # Now locate outmode
            found = False
            for x in range(0, len(elems)):
                if elems[x] == outmode:
                    found = True
                    new = elems[x + 1]
                    break
            if not found:
                continue
            conv[key][old] = new
    return conv


def _checkpdb(pdb_in, pdb_out, conversion):
    """
    Checks if all atom names are already in the PDB file

    Parameters
    ----------
    pdb_in : PDBFile
      the pdb file
    pdb_out : PDBFile
      the output pdb file
    conversion :
      a dictionary with conversion

    Returns
    -------
    bool
      True if all atom names are already in the PDB file
      False if not
    """

    for resnum in pdb_in.residues:
        residue = pdb_in.residues[resnum]

        if residue.name.upper(
        ) == "HID":  # Valid only when the pdb input file is from AMBER.
            resname = "HIS"
            pdb_out.residues[resnum].name = resname
        else:
            resname = residue.name.upper()

        for atomnum in range(len(residue.atoms)):
            atomname = residue.atoms[atomnum].name.upper()
            # Check if the atom name is either in the backbone
            # or a sidechain residue
            if not (atomname in conversion["backbone"].values()
                    or atomname in conversion[resname].values()):
                return False
    return True


def pdb2pms(pdb_in, forcefield, conversion_file):
    """
    Convert atom names to ProtoMS standard

    The protein object is not modified by this routine, but the Residue
    and Atom objects are.

    Parameters
    ----------
    pdb_in : PDBFile
      the pdb file to modify inline
    forcefield : string
      the style of the input PDB-file
    conversion_file :
      a file with conversion rules

    Returns
    -------
    a PDBFile instance with ProtoMS names
    """

    logger.debug("Running pdb2pms with arguments: ")
    logger.debug("\tpdb_in          = %s" % pdb_in)
    logger.debug("\tforcefield      = %s" % forcefield)
    logger.debug("\tconversion_file = %s" % conversion_file)
    logger.debug("This will rename atoms in a PDB-file to match "
                 "ProtoMS naming convention")

    pdb_out = pdb_in.copy()
    conversion = read_convfile(
        conversion_file, inmode=forcefield, outmode="protoms")
    if _checkpdb(pdb_in, pdb_out, conversion):
        logger.info(
            "It seems that %s already contains the correct ProtoMS names."
            " Aborting the conversion." % pdb_in.name)
        return pdb_in
    for resnum in pdb_in.residues:
        residue = pdb_in.residues[resnum]
        if residue.name.upper(
        ) == "HID":  # Valid only when the pdb input file is from AMBER.
            resname = "HIS"
            pdb_out.residues[resnum].name = resname
        else:
            resname = residue.name.upper()
        for atomnum in range(len(residue.atoms)):
            atomname = residue.atoms[atomnum].name.upper()
            if atomname in conversion["backbone"]:
                newName = conversion["backbone"][atomname]
            elif atomname in conversion[resname]:
                newName = conversion[resname][atomname]
            elif atomname[1:len(atomname)] + atomname[0] in conversion[
                    resname]:
                newName = conversion[resname][(
                    atomname[1:len(atomname)] + atomname[0])]
            elif atomname[-1] + atomname[0:(len(atomname) - 1)] in conversion[
                    resname]:
                newName = conversion[resname][(
                    atomname[-1] + atomname[0:(len(atomname) - 1)])]
            else:
                newName = atomname
                logger.warning(
                    "Warning: atom %s in residue %s %i not found in %s. This "
                    "atom has not been converted. Check atom names are validb."
                    % (atomname, resname, resnum, conversion_file))
            pdb_out.residues[resnum].atoms[atomnum].name = newName
            pdb_out.residues[resnum].atoms[atomnum].resname = resname
    return pdb_out


def make_dummy(pdbfile):
    """
    Make a dummy PDB structure to match a molecule

    It will place the dummy at the centre of the molecule

    Parameters
    ----------
    pdbfile : string or PDBFile
      the pdb structure of the molecule

    Returns
    -------
    PDBFile instance
      the dummy structure
    """

    logger.debug("Running make_dummy with arguments: ")
    logger.debug("\tpdbfile = %s" % pdbfile)
    logger.debug("This will make a dummy for the molecule in pdbfile")

    if isinstance(pdbfile, six.string_types):
        pdbfile = simulationobjects.PDBFile(filename=pdbfile)

    # HEADER dummy
    # ATOM      1  C1  DDD     1       0.000   0.000   0.000  1.00  0.00

    atom = simulationobjects.Atom()
    atom.index = 1
    atom.resindex = 1
    atom.name = " C1 "
    atom.resname = "DDD"
    atom.coords = np.array(pdbfile.getCenter(), copy=True)

    residue = simulationobjects.Residue()
    residue.index = 1
    residue.name = "DDD"
    residue.atoms.append(atom)

    dummy = simulationobjects.PDBFile()
    dummy.residues[1] = residue
    dummy.header = "HEADER dummy\n"

    return dummy
