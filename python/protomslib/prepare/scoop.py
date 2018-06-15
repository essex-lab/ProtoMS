import itertools
import logging
import operator
import os
import six
import numpy as np
from .. import simulationobjects

logger = logging.getLogger('protoms')


def scoop(protein,
          ligand,
          innercut=16,
          outercut=20,
          flexin='full',
          flexout='sidechain',
          terminal='neutralize',
          excluded=[],
          added=[],
          scooplimit=10):
    """
    Generates a scoop from protein structure

    The protein object is not modified by this routine.

    Parameters
    ----------
    protein : PDBFile
        pdb instance for the protein to be scooped
    ligand : PDBFile or string
        Either pdb instance for ligand to be scooped around, or string giving
        the name of a file containing 3 floats to act as coords for scoop
        centre, or a string with three floating point numbers
    innercut : float, optional
        Maximum distance from ligand defining inner region of the scoop
    outercut : float, optional
        Maximum distance from ligand defining outer region of the scoop
    flexin : string, optional
        Gives the degree of flexibility for residues of the inner region
        Can be 'rigid', 'sidechain' or 'flexible'
    flexout : string, optional
        As flexin but for residues of the outer scoop
    terminal : string, optional
        Determines what to do with terminal residues
        Can be 'keep', 'doublekeep','neutralize'
    excluded : list of int
        List of indices for residues to be excluded from scoop
    added : list of int
        List of indices for residues to be included in outer scoop

    Returns
    -------
    PDBFile
        an object representing the scooped protein
    """

    logger.debug("Running scoop with arguments: ")
    logger.debug("\tprotein  = %s" % protein)
    logger.debug("\tligand   = %s" % ligand)
    logger.debug("\tinnercut = %f" % innercut)
    logger.debug("\toutercut = %f" % outercut)
    logger.debug("\tflexin   = %s" % flexin)
    logger.debug("\tflexout  = %s" % flexout)
    logger.debug("\tterminal  = %s" % terminal)
    logger.debug("\texcluded = %s" % " ".join("%d" % e for e in excluded))
    logger.debug("\tadded    = %s" % " ".join("%d" % a for a in added))
    logger.debug("This will generate a truncated version for a protein")

    pdb_out = protein.copy()

    centerarray = None
    if isinstance(ligand, six.string_types):
        if os.path.isfile(ligand):
            centerarray = np.loadtxt(ligand)
        else:
            centerarray = np.array(ligand.strip().split(), dtype=float)

    # Either if a string or a file with center coordinates was passed,
    # we will make it into a dummy pdb object
    if centerarray is not None:
        ligand = simulationobjects.PDBFile()
        ligand.residues = {0: simulationobjects.Residue()}
        ligand.residues[0].addAtom(simulationobjects.Atom(coords=centerarray))

    assert flexin in ['sidechain', 'flexible', 'rigid']
    assert flexout in ['sidechain', 'flexible', 'rigid']

    if len(ligand.residues) > 1:
        print("More than one ligand in input. Scooping around everything...")

    # Build inner and outer lists
    # If any heavy atom of a residue falls within the cutoff distance of
    # any atom of the ligand, add then to the appropriate list
    outer = []
    inner = []
    kill_list = [True] * len(pdb_out.residues)
    in_kill_list = [True] * len(pdb_out.residues)
    for res in pdb_out.residues:
        for atom in pdb_out.residues[res].atoms:
            if (atom.name.startswith(('H', 'h', '1'))
                    or atom.name in ('N', 'C', 'O')):
                continue
            for lig in six.itervalues(ligand.residues):
                for lat in lig.atoms:
                    if lat.name.startswith(('H', 'h')):
                        continue
                    distance = np.linalg.norm(atom.coords - lat.coords)
                    if distance < outercut:
                        kill_list[res - 1] = False
                    if distance < innercut:
                        in_kill_list[res - 1] = False
        if not kill_list[res - 1] and \
           in_kill_list[res - 1] and res not in excluded:
            outer.append(res)
        if kill_list[res - 1] and res in added:
            outer.append(res)
        if not in_kill_list[res - 1]:
            inner.append(res)

    kill_list = np.array(kill_list)
    in_kill_list = np.array(in_kill_list)

    # Also eliminate Xray waters beyond outer cutoff
    waters = []
    for mol in pdb_out.solvents:
        kill = True
        for atom in pdb_out.solvents[mol].atoms:
            for lig in six.itervalues(ligand.residues):
                for lat in lig.atoms:
                    distance = np.linalg.norm(atom.coords - lat.coords)
                    if distance < outercut:
                        kill = False
                        break
        if not kill:
            waters.append(mol)

    both = sorted(inner + outer)

    # All CYX residues and the backbone of neighbouring residues
    # must be fixed to preserve disulphide bridges.
    rigid = []
    backBoneRigid = []
    for res in both:
        if pdb_out.residues[res].name == 'CYX':
            rigid.append(res)
            if res != 0:
                backBoneRigid.append(res - 1)
            if res != len(pdb_out.residues):
                backBoneRigid.append(res + 1)

    # Constrain outer residues only if the number of discarded
    # residues is less than scooplimit.
    if kill_list.sum() > scooplimit:
        if flexout in ['rigid', 'sidechain']:
            backBoneRigid += [res for res in outer]
        if flexout == 'rigid':
            rigid += [res for res in outer]
    # Same thing for inner residues
        if flexin in ['rigid', 'sidechain']:
            backBoneRigid += [res for res in inner]
        if flexin == 'rigid':
            rigid += [res for res in inner]

    outres = []
    rigidBB = []
    rigidSC = []
    count = 0
    for res in both:
        outres.append(pdb_out.residues[res])
        count += 1
        if res in backBoneRigid:
            rigidBB.append(count)
        if res in rigid:
            rigidSC.append(count)

    # Need to turn residue lists into nice string of residue ranges for
    # ProtoMS to interpret. Not crazy about this but works.
    outBB = ''
    for k, g in itertools.groupby(
            enumerate(rigidBB), key=lambda x: x[0] - x[1]):
        r = list(map(operator.itemgetter(1), g))
        if len(r) > 1:
            outBB += '%d-%d, ' % (min(r), max(r))
        else:
            outBB += '%d, ' % r[0]
    outBB = outBB[:-2]

    outSC = ''
    for k, g in itertools.groupby(
            enumerate(rigidSC), key=lambda x: x[0] - x[1]):
        r = map(operator.itemgetter(1), g)
        if len(r) > 1:
            outSC += '%d-%d, ' % (min(r), max(r))
        else:
            outSC += '%d, ' % r[0]
    outSC = outSC[:-2]

    # Purge residues outside the outer scoop from the protein pdb and save it
    # Only do it if the number of discarded residues is less than scooplimit.
    if kill_list.sum() > scooplimit:
        for res in list(pdb_out.residues.keys()):
            if res not in both:
                pdb_out.residues.pop(res)
        for res in list(pdb_out.solvents.keys()):
            if res not in waters:
                pdb_out.solvents.pop(res)
    else:
        # print "Not scooping, lower than limit."
        logger.info(
            "Not scooping. Number of residues removed from the protein is"
            " too small (%s less than %s)" % (kill_list.sum(), scooplimit))

    # Check of terminal residue
    headerscoop = False
    if terminal != 'keep':
        nres = pdb_out.residues[sorted(pdb_out.residues.keys())[0]]
        cres = pdb_out.residues[sorted(pdb_out.residues.keys())[-1]]

        # Check if the N terminal is charged,
        # only works for now with ProtoMS names
        charged_nres = False
        for atom in nres.atoms:
            if atom.name.upper().strip() in ["HN1", "HN2", "HN3"]:
                charged_nres = True

        # Check if the C terminal is charged
        charged_cres = False
        for atom in cres.atoms:
            if atom.name.upper().strip() in ["OT", "OXT"]:
                charged_cres = True

        # If both are charged and terminal is doublekeep,
        # we will not tell ProtoMS this is a scoop and keep them
        if charged_nres and charged_cres and terminal == 'doublekeep':
            headerscoop = False
        # Otherwise we will neutralize and tell ProtoMS this is a scoop
        else:
            headerscoop = True
            if charged_nres:
                remthis = []
                for atom in nres.atoms:
                    if atom.name.upper().strip() in ["HN1", "HN2", "HN3"]:
                        remthis.append(atom)
                for atom in remthis:
                    nres.atoms.remove(atom)
            if charged_cres:
                remthis = []
                for atom in cres.atoms:
                    if atom.name.upper().strip() in ["OT", "OXT"]:
                        remthis.append(atom)
                for atom in remthis:
                    cres.atoms.remove(atom)

    header = "REMARK Original file: %s\n" % pdb_out
    if kill_list.sum() > scooplimit:
        header += "REMARK Scoop of %s was made\n" % pdb_out
        header += "REMARK Inner Region : %8.2f Angstrom radius\n" % innercut
        header += "REMARK Outer Region : %8.2f Angstrom radius\n" % outercut
        header += "REMARK Num Residues %d ( %d inner %d outer )\n" % (
            len(inner) + len(outer), len(inner), len(outer))
        header += "REMARK %d residues have a fixed backbone\n" % (len(rigidBB))
        header += "REMARK %d residues are fixed\n" % (len(rigidSC))
        header += "REMARK flexibility of the inner part : %s\n" % flexin
        header += "REMARK flexibility of the outer part : %s\n" % flexout
        header += "REMARK ProtoMS keyword to use\n"


# check that the outBB and outSC are short enough to be printed on one line
    linelength = 270
    if len(outBB) <= linelength:
        if not len(rigidBB) == 0:
            header += "REMARK chunk fixbackbone 1 %s\n" % outBB
    else:
        print("REMARK chunk fixbackbone line is too long for ProtoMS, "
              "please split over two lines")
    if len(outSC) <= linelength:
        if not len(rigidSC) == 0:
            header += "REMARK chunk fixresidues 1 %s\n" % outSC
    else:
        print("REMARK chunk fixresidues line is too long for ProtoMS, "
              "please split over two lines")

    if kill_list.sum() > scooplimit:
        header += "REMARK Xray Water within %8.2f Angstrom\n" % outercut
        header += "REMARK of the ligand\n"
        if headerscoop:
            header += 'HEADER scoop\n'

    pdb_out.header = header

    return pdb_out
