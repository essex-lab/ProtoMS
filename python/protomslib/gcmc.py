from __future__ import print_function
import logging
import numpy as np
import six
from . import simulationobjects
from .prepare.water import convertwater, rotatesolute

logger = logging.getLogger("protoms")


def make_gcmcbox(pdb, filename, padding=2.0):
    """
    Make a GCMC/JAWS-1 simulation box around a PDB-structure

    Parameters
    ----------
    pdb : PDBFile object
      the PDB structure
    filename : string
      the name of the output file
    padding : float, optional
      the amount of extra space around the ligand to add
    """

    logger.debug("Running make_gcmcbox with arguments: ")
    logger.debug("\tpdb      = %s" % pdb)
    logger.debug("\tfilename = %s" % filename)
    logger.debug("\tpadding  = %d" % padding)
    logger.debug("This will make a simulation box for GCMC/JAWS-1")

    # Create a box around the solute and pad it with two Angstromgs
    box = pdb.getBox()
    box["origin"] = box["origin"] - padding
    box["len"] = box["len"] + 2.0 * padding

    return box


def print_bequil(boxlen):
    betamu = -10.47
    boxvolume = boxlen[0] * boxlen[1] * boxlen[2]
    bequil = betamu + np.log(boxvolume / 30.0)
    print("Volume of GCMC box:", np.round(boxvolume, 2))
    print("Bequil:", np.round(bequil, 2))


def _create_res(
    atnames=["O00"],
    resname="sol",
    positions=[np.array([0.0, 0.0, 0.0])],
    resind=1,
):
    """
    Create a residue object

    Parameters
    ----------
    atnames : list of strings
      a list containing the names of all
    the atoms in the residue
      resname : string
    the name of the residue
      positions : list of numpy arrays
    a list with the positions of all
      the atoms in the residue

    Returns
    -------
    Residue object
      the residue built with the given atoms
    """
    resobj = simulationobjects.Residue()
    resobj.name = resname
    resobj.index = resind
    for ind, atom in enumerate(atnames):
        atobj = simulationobjects.Atom()
        atobj.name = atom
        atobj.resname = resname
        atobj.index = ind
        atobj.resindex = resind
        try:
            atobj.coords = positions[ind]
        except Exception:
            atobj.coords = np.array([0.0, 0.0, 0.0])
        resobj.addAtom(atobj)
    return resobj


def distribute_particles(
    box, particles, watermodel="t4p", resname="WA1", partnumb=None
):
    """
    Randomly distribute molecules in a box

    Parameters
    ----------
    box : a list or dictionary
      the origin and length of the box
      where the molecules will be placed
    particles : string or PDB object
      the number of waters to include in the box
      or the pdb object or filename with the molecules
      to include in the box
    watermodel : string, optional
      (only used when particles is a number)
      for instance, "t4p", "tip4p", "t3p" "tip3p",
      the water model for the generated waters
    partnumb : string,optional
      (only used when particles is a file)
      the number of particles. If not specified it is set
      to be as many as there are in the file

    Returns
    -------
    pdb object
      the pdb object with the molecules
      distributed randomly in the box
    """

    logger.debug("Running distribute_particles with arguments: ")
    if isinstance(box, dict):
        tmpstr = "\tbox = "
        for key in box:
            tmpstr = tmpstr + "%s %s, " % (key, str(box[key]))
        logger.debug(tmpstr.strip(" ,"))
    else:
        logger.debug("\tbox = %s" % " ".join(box))
    if isinstance(particles, str):
        logger.debug("\tparticles = %s" % particles)
    else:
        logger.debug("\tparticles %s" % particles.name)
    logger.debug("\twatermodel %s" % watermodel)
    logger.debug("\tpartnumb %s" % partnumb)
    logger.debug(
        "This will distribute the molecules given in 'particles' "
        "randomly in the box given in box. If a number is given in"
        " 'particles', it assumes them to be waters."
    )

    if isinstance(box, list):
        try:
            box = [float(val) for val in box]
        except Exception:
            raise simulationobjects.SetupError(
                "The box dimensions %s could not be correctly interpreted"
                % box
            )
        orig = np.array(box[:3])
        length = np.array(box[3:])
    elif isinstance(box, dict):
        orig = box["origin"]
        length = box["len"]
    else:
        # Ensure error raised for unknown type handling
        raise TypeError(box)

    if isinstance(particles, str) and particles.isdigit():
        watnumb = int(particles)
        particles = simulationobjects.PDBFile()
        for i in range(1, watnumb + 1):
            pad = 0.3  # Water oxygens are at least this far away (in Angs.) from the edge of the box
            oxpos = np.random.rand(len(length)) * np.array(
                length - 2 * pad
            ) + np.array(orig + pad)
            particles.solvents[i] = _create_res(
                resname=resname, positions=[oxpos]
            )
        particles = convertwater(
            particles, watermodel, "y", watresname=resname
        )

    else:
        if isinstance(particles, str):
            try:
                pdbobj = simulationobjects.PDBFile()
                pdbobj.read(particles)
            except Exception:
                raise simulationobjects.SetupError(
                    "The pdb file %s could not be found" % particles
                )
            particles = pdbobj
        if particles.residues:
            parts = particles.residues
        elif particles.solvents:
            parts = particles.solvents
        else:
            raise simulationobjects.SetupError(
                "No molecule could be found in %s" % particles
            )
        if partnumb is not None:
            particles = simulationobjects.PDBFile()
            for i in range(1, int(partnumb) + 1):
                allcoords = [
                    myat.coords for myat in parts[parts.keys()[0]].atoms
                ]
                allnames = [myat.name for myat in parts[parts.keys()[0]].atoms]
                particles.residues[i] = _create_res(
                    resname=parts[parts.keys()[0]].name,
                    positions=allcoords,
                    atnames=allnames,
                    resind=i,
                )
            parts = particles.residues
        for ind, keypart in enumerate(parts):
            parts[keypart].getCenter()
            displace = [
                coord_len * np.random.uniform()
                + orig[jnd]
                - parts[keypart].center[jnd]
                for jnd, coord_len in enumerate(length)
            ]
            rotated_coords = rotatesolute(
                np.array([myat.coords for myat in parts[keypart].atoms]),
                np.random.uniform(0, 2 * np.pi),
                np.random.uniform(0, 2 * np.pi),
                np.random.uniform(0, 2 * np.pi),
            )
            for jnd, atom in enumerate(parts[keypart].atoms):
                atom.coords = np.array(
                    [
                        rotated_coords.item((jnd, i)) + displace[i]
                        for i in range(3)
                    ]
                )

    for ind, coord in enumerate(orig):
        h_parts = particles.header.strip().split()
        particles.header = "%s %.4f %s %.4f " % (
            " ".join(h_parts[:ind]),
            coord,
            " ".join(h_parts[ind : ind * 2]),
            coord + length[ind],
        )
    particles.header = "HEADER box %s\n" % particles.header
    return particles


def clear_gcmcbox(gcmcbox, waters):
    """
    Removes solvent molecules from the GCMC/JAWS-1 box as they will not be
    able to move during the simulation

    Parameters
    ----------
    gcmcbox : string or PDBFile object
      the gcmcbox
    waters : string or PDBFile object
      the water molecule to check

    Returns
    -------
    int
      the number of removed water molecules
    PDBFile
      the cleared solvation box
    """

    logger.debug("Running clear_gcmcbox with arguments: ")
    logger.debug("\tgcmcbox = %s" % gcmcbox)
    logger.debug("\twaters = %s" % waters)
    logger.debug("This will remove solvent molecules within the GCMC/JAWS box")

    # The amount to clear in excess of the box limits, as only the
    # oxygen atoms of the water molecules are used to decide
    # whether the water molecule is in the box or not.
    extend = 1.0
    if isinstance(gcmcbox, six.string_types):
        gcmcbox = simulationobjects.PDBFile(filename=gcmcbox)
    if isinstance(waters, six.string_types):
        waters = simulationobjects.PDBFile(filename=waters)

    # Try to find box information in the header
    box = simulationobjects.find_box(gcmcbox)
    if box is None:
        # if that fails, take the extent of the PDB structure
        box = gcmcbox.getBox()
    if "origin" not in box:
        box["origin"] = box["center"] - box["len"] / 2.0
    box_max = box["origin"] + box["len"]
    box_min = box["origin"]

    # Remove waters that are inside the GCMC/JAWS-1 box
    nrem = 0
    removethese = []
    for soli in waters.solvents:
        xyz = waters.solvents[soli].atoms[0].coords
        if np.all(xyz < (box_max + extend)) and np.all(
            xyz > (box_min - extend)
        ):
            logger.debug("Removing water %d from %s" % (soli, waters))
            nrem = nrem + 1
            removethese.append(soli)
    for soli in removethese:
        del waters.solvents[soli]
    logger.info(
        "Removed %d water molecules from %s that were inside "
        "the GCMC/JAWS box %s" % (nrem, waters, gcmcbox)
    )

    return nrem, waters
