# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Routines to make ProtoMS single topology template

This module defines these public function
make_single
write_map
summarize

Can be executed from the command line as a stand-alone program
"""

import copy
import logging

import numpy as np

import simulationobjects as sim

logger = logging.getLogger('protoms')


def _make_dict(atoms, moltem, pdbres, objdict=None, onlypdb=False):
    """ 
    Finds a TemplateAtom and an Atom object for each atom in atoms

    Parameters
    ---------- 
    atoms : list of string
      the atom to lookup
    moltem : MolTemplate
      the template to search in
    pdbres : Residue
      the pdb residue to search in
    objdict : dictionary, optional
      dictionary to add information to
    onlypdb : boolean, optional
      flag indicating if only find Atom object

    Returns
    -------
    dictionary
      for each atom in atoms parameter, a dictionary with
      keys "tem" and "pdb", which has a TemplateAtom and Atom object
      as its value
    """
    if objdict is None:
        objdict = {}
    for atom in atoms:
        if atom in objdict:
            continue
        objdict[atom] = {}
        if not onlypdb:
            for tatom in moltem.atoms:
                if tatom.name.strip().upper() == atom:
                    objdict[atom]["tem"] = tatom
                    break
        for patom in pdbres.atoms:
            if patom.name.strip().upper() == atom:
                objdict[atom]["pdb"] = patom
                break
        if "pdb" not in objdict[atom]:
            logger.warning(
                "Warning: Could not find atom %s in the pdb-file, cannot look for map for this atom" % atom)
            objdict[atom]["pdb"] = None
    return objdict


def _make_map(tem1, tem2, pdb1, pdb2, cmap):
    """
    Make a correspondance map from template 1 to template 2

    The cmap dictionary will be populated with key value pairs:
    cmap[atom1] = atom2,
    where atom1 belongs to tem1 and atom2 belongs to tem2
    atom2 can be dum, in case non-correspondance is assumed

    Parameters
    ----------
    tem1 : TemplateFile
      first template file
    tem2 : TemplateFile
      second template file
    pdb1 : PDBFile
      structure template for tem1
    pdb2 : PDBFile
      structure template for tem2
    cmap : dictionary of string
      full, partial or empty map of correspondance

    Raises
    ------
    SetupError
      if dummies in cmap.keys()
    """

    # Create a list of atom names in both templates
    # The purpose of the routine is then to empty these lists!
    not_taken1 = [atom.name.strip().upper()
                  for atom in tem1.templates[0].atoms]
    not_taken2 = [atom.name.strip().upper()
                  for atom in tem2.templates[0].atoms]

    # First we will look in the given cmap to see if we already have some maps
    # given
    for atom1 in cmap:
        atom2 = cmap[atom1]
        if atom1.startswith("DUM"):
            sim.SetupError("Do not support dummies in V0 at the moment.")
        if atom1 in not_taken1:
            if atom2 != "DUM":
                if atom2 in not_taken2:
                    not_taken1.remove(atom1)
                    not_taken2.remove(atom2)
                else:
                    logger.warning(
                        "Warning: Could not find atom %s in the second template file, will ignore this map." % atom2)
                    del cmap[atom1]
            else:
                not_taken1.remove(atom1)
        else:
            logger.warning(
                "Warning: Could not find atom %s in the first template file, will ignore this map." % atom2)
            del cmap[atom1]
    logger.info("Pre-defined maps:\n%s" % " ".join("%s-%s" %
                                                   (atom1, cmap[atom1]) for atom1 in sorted(cmap.keys())))

    # Next we will calculate pair-wise distance from the given pdb structures,
    # use a cut-off of 0.02 A^2 to determine possible pair and then use atom
    # type as the final filter
    if not_taken1 or not_taken2:
        objdict1 = _make_dict(not_taken1, tem1.templates[0], pdb1.residues[1])
        objdict2 = _make_dict(not_taken2, tem2.templates[0], pdb2.residues[1])
        found = {}
        for atom1 in not_taken1:
            if objdict1[atom1]["pdb"] is None:
                continue
            for atom2 in not_taken2:
                if objdict2[atom2]["pdb"] is None:
                    continue
                dist2 = objdict1[atom1]["pdb"].coords - \
                    objdict2[atom2]["pdb"].coords
                dist2 = np.sum(dist2 * dist2)
                # The atom type is the 0:th parameter of a clj parameter
                type1 = objdict1[atom1]["tem"].param0.params[0]
                type2 = objdict2[atom2]["tem"].param0.params[0]
                if dist2 < 0.02 and type1 == type2:
                    found[atom1] = atom2
                    break
        logger.info("Distance/atom type maps:\n%s" % " ".join("%s-%s" %
                                                              (atom1, found[atom1]) for atom1 in sorted(found.keys())))
        for atom1 in found:
            atom2 = found[atom1]
            not_taken1.remove(atom1)
            not_taken2.remove(atom2)
            cmap[atom1] = atom2

    # Now we have gotten so-far that we have to ask the user for the rest...

    if not_taken1 or not_taken2:

        # Print out useful information
        logger.info("")
        logger.info("These are the un-matched atoms of template 1: %s" %
                    " ".join(not_taken1))
        logger.info("These are the un-matched atoms of template 2: %s" %
                    " ".join(not_taken2))

        logger.info("")
        logger.info("These are the distances (A): ")
        logger.info("%8s%s" % ("", "".join("%8s" %
                                           atom for atom in not_taken1)))
        for atom2 in not_taken2:
            outstr = "%8s" % atom2
            for atom1 in not_taken1:
                dist = objdict1[atom1]["pdb"].coords - \
                    objdict2[atom2]["pdb"].coords
                dist = np.sqrt(np.sum(dist * dist))
                outstr = outstr + "%8.3f" % dist
            logger.info(outstr)
        logger.info("")

        found = {}
        for atom1 in not_taken1:
            atom2 = ""
            dummy = False
            while atom2 not in not_taken2 and not dummy:
                print "Enter the corresponding atom for %s: " % atom1,
                test = raw_input().strip().upper()
                if len(test) == 0 or test.startswith("DUM"):
                    dummy = True
                    atom2 = ""
                else:
                    atom2 = test
            if dummy:
                found[atom1] = "DUM"
            else:
                found[atom1] = atom2
        logger.info("")
        logger.info("User maps:\n%s" % " ".join("%s-%s" %
                                                (atom1, found[atom1]) for atom1 in sorted(found.keys())))
        for atom1 in found:
            atom2 = found[atom1]
            not_taken1.remove(atom1)
            if atom2 != "DUM":
                not_taken2.remove(atom2)
            cmap[atom1] = atom2


def _make_ele_tem(tem1, tem2, cmap):
    """
    Make new single topology template file for electrostatic leg

    Parameters
    ----------
    tem1 : TemplateFile
      first template file
    tem2 : TemplateFile
      second template file
    cmap : dictionary of string
      map of correspondance

    Returns
    -------
    TemplateFile
      the created template file
    """

    # Make copies of the objects first so we can manipulate them without
    # modifying the original templates
    tem1 = copy.deepcopy(tem1)
    tem2 = copy.deepcopy(tem2)

    # Copy all clj parameters from tem2 to tem1
    if tem1.cljparams:
        start = max(p.index for p in tem1.cljparams) + 1
    else:
        start = -50000
    for index, param in enumerate(tem2.cljparams, start):
        if index > 0:
            param.index = index
        tem1.cljparams.append(param)

    # Change the parm1 for all atoms in tem1
    nextidx = max(p.index for p in tem1.cljparams) + 1
    for atom1 in tem1.templates[0].atoms:
        atom2 = cmap[atom1.name.upper()]
        if atom2 != "DUM":  # Needs to find atom in tem2
            for atom3 in tem2.templates[0].atoms:
                if atom3.name.strip().upper() == atom2:
                    atom1.param1 = atom3.param0
                    tempch = atom1.param1.params[2]
                    # This make sure that we do not change vdw
                    atom1.param1.params = list(atom1.param0.params)
                    atom1.param1.params[2] = tempch
                    break
        else:
            # Needs to create a new parameter
            newparam = copy.deepcopy(atom1.param1)
            newparam.index = nextidx
            nextidx = nextidx + 1
            # For clj parameters, 3:rd parameter is the charge
            newparam.params[2] = "0.00000"
            atom1.param1 = newparam
            tem1.cljparams.append(newparam)

    return tem1


def _make_vdw_tem(tem1, tem2, pdb1, pdb2, cmap, usepdb=True, gaffversion="gaff16"):
    """
    Make new single topology template file for van der Waals leg

    Parameters
    ----------
    tem1 : TemplateFile
      first template file
    tem2 : TemplateFile
      second template file
    pdb1 : PDBFile
      structure template for tem1
    pdb2 : PDBFile
      structure template for tem2
    cmap : dictionary of string
      map of correspondance
    usepdb : boolean, optional
      flag indicating if using variable geometries from pdb-file
    gaffversion : str, optional
      version of gaff to use, defaults to gaff16

    Returns
    -------
    TemplateFile
      the created template file
    """

    def find_param(atoms, temset, gaffset):
        """ Find force field index and equilibrium value for atom types
        """
        ratoms = atoms[::-1]
        idx = -1
        # Try to look in the temset first
        for atoms2 in temset:
            if atoms2.atoms == atoms or atoms2 == ratoms:
                return atoms2.param.index, atom2.param.params[1]
        # Next loop through the GAFF set
        if len(atoms) == 4:
            param = gaffset.get_params(atoms)
            return param.index, param.terms
        else:
            param = gaffset.get_params(atoms)
            return param.index, param.b0
        return -1, -1

    def find_pdbparam(atoms, objdict):
        """ Find bond lengths and angles in pdb files
        """
        allcoords = [objdict[atom]["pdb"].coords for atom in atoms]
        allcoords = np.array(allcoords)
        if len(atoms) == 2:
            return np.sqrt(np.sum((allcoords[0, :] - allcoords[1, :])**2))
        elif len(atoms) == 3:
            return sim.angle_atms(allcoords[0, :], allcoords[1, :], allcoords[2, :]) * 180.0 / np.pi
        return 0.0

    def make_variable(atom):
        """ Add variable geometry to a molecular template
        """
        # Skip this if the atom is bonded to a dummy atom
        if atom.bondedto in ["DM1", "DM2", "DM3"]:
            return
        # This flags the creating of a dummy angle
        defangle = not atom.angleto in ["DM1", "DM2", "DM3"]

        # Store away z-matrix atoms and their names
        zmatatoms = [atom, atom.bondedto]
        if defangle:
            zmatatoms.append(atom.angleto)
        zmatnames1 = [zatom.name.upper() for zatom in zmatatoms]
        zmatnames2 = [cmap[zatom.name.upper()] for zatom in zmatatoms]
        atomtypes1 = [zatom.param0.params[0]
                      for zatom in zmatatoms]  # Atom types for template 1
        if usepdb:  # If we should take the geometry from the pdb-files
            _make_dict(zmatnames1, tem1.templates[0], pdb1.residues[
                       1], objdict1, onlypdb=True)  # Find the Atom object and put in a dictionary
            bond1 = find_pdbparam(zmatnames1[:2], objdict1)
        else:  # If we should take the geometry from GAFF
            dummy, bond1 = find_param(atomtypes1[:2], temsets[
                                      "bond"], gaffsets["bond"])
        if isinstance(atom.param1, int):  # Check if this atom should be perturbed to a dummy
            tem1.templates[0].variables.append(
                "# %s to dummy" % "-".join(atomtypes0[:2]))
            tem1.templates[0].variables.append("variable %s %s bond %.3f %.3f" % (
                atom.name, atom.residue, bond1, 0.200))  # Shrink bond to within vdw-sphere
        else:  # State V1 has cljparameters
            atomtypes2 = [zatom.param1.params[0] for zatom in zmatatoms]
            if usepdb:
                _make_dict(zmatnames2, tem2.templates[
                           0], pdb2.residues[1], objdict2, onlypdb=True)
                bond2 = find_pdbparam(zmatnames2[:2], objdict2)
                if defangle:
                    angle1 = find_pdbparam(zmatnames1, objdict1)
                    angle2 = find_pdbparam(zmatnames2, objdict2)
            else:
                dummy, bond2 = find_param(atomtypes2[:2], temsets[
                                          "bond"], gaffsets["bond"])
                if defangle:
                    dummy, angle1 = find_param(
                        atomtypes1, temsets["angle"], gaffsets["angle"])
                    dummy, angle2 = find_param(
                        atomtypes2, temsets["angle"], gaffsets["angle"])
            # Add variable geometry for bond and angle
            if bond1 != bond2:
                tem1.templates[0].variables.append("# %s to %s at atoms %s" % (
                    "-".join(atomtypes1[:2]), "-".join(atomtypes2[:2]), "-".join(zmatnames1[:2])))
                tem1.templates[0].variables.append(
                    "variable %s %s bond %.3f %.3f" % (atom.name, atom.residue, bond1, bond2))
            if defangle and angle1 != angle2:
                tem1.templates[0].variables.append("# %s to %s at atoms %s" % (
                    "-".join(atomtypes1), "-".join(atomtypes2), "-".join(zmatnames1)))
                tem1.templates[0].variables.append(
                    "variable %s %s angle %.3f %.3f" % (atom.name, atom.residue, angle1, angle2))

    # MAIN routine

    # Make copies of the objects first so we can manipulate them without
    # modifying the original templates
    tem1 = copy.deepcopy(tem1)
    tem2 = copy.deepcopy(tem2)

    # Copy all clj parameters from tem2 to tem1, both of the sets will be
    # modified below
    if tem1.cljparams:
        start = max(p.index for p in tem1.cljparams) + 1
    else:
        start = -50000
    for index, param in enumerate(tem2.cljparams, start):
        if index > 0:
            param.index = index
        tem1.cljparams.append(param)

    # Change the parm0 and parm1 for all atoms in tem1
    atomlist = []  # Will contain atoms where the atomtype changes
    nextidx = max(p.index for p in tem1.cljparams) + 1
    for atom1 in tem1.templates[0].atoms:
        atom2 = cmap[atom1.name.upper()]
        if atom2 != "DUM":  # Needs to find atom in tem2
            for atom3 in tem2.templates[0].atoms:
                if atom3.name.strip().upper() == atom2:
                    atom1.param0.params[2] = atom3.param0.params[
                        2]  # Change the charges
                    atom1.param1 = atom3.param0  # Change the param1, charge AND vdw
                    # If the atom type is different, store away it
                    if atom1.param0.params[0] != atom1.param1.params[0]:
                        atomlist.append(atom1)
                    break
        else:
            # Needs to create a new parameter
            newparam = copy.deepcopy(atom1.param1)
            newparam.index = nextidx
            nextidx = nextidx + 1
            tem1.cljparams.append(newparam)
            newparam.params[2] = "0.00000"  # Set the charge to zero
            atom1.param0 = newparam
            atom1.param1 = 100  # This is a dummy parameter
            atomlist.append(atom1)

    # Now we need to change the connectivity of the template so that the intramolecular energies/geometry are correctly calculated
    # print "\nThe following atoms change type: %s"%" ".join("%s"%atom.name
    # for atom in atomlist)

    # Setup parameter sets to search
    gaffname = sim.standard_filename(gaffversion + ".ff", "parameter")
    gaffsets = {ptype: sim.ParameterSet(ptype, gaffname) for ptype in [
        "bond", "angle", "dihedral"]}
    temsets = {"bond": tem1.bondatoms, "angle": tem1.angleatoms,
               "dihedral": tem1.dihedralatoms}

    # Now loop over all atoms that changes type
    objdict1 = _make_dict([atom.name for atom in atomlist], tem1.templates[
                          0], pdb1.residues[1], onlypdb=True)
    objdict2 = _make_dict([cmap[atom.name.upper()] for atom in atomlist if cmap[
                          atom.name.upper()] != "DUM"], tem2.templates[0], pdb2.residues[1], onlypdb=True)
    variablemade = [False] * len(tem1.templates[0].atoms)
    newdihedrals = []
    for atom in atomlist:
        # Add variable geometry to this atom
        make_variable(atom)
        # Then add variable geometry for all other atoms bonded to this atom,
        # except if it is in the atomlist
        for i, atom2 in enumerate(tem1.templates[0].atoms):
            if variablemade[i] or atom2 in atomlist:
                continue
            if atom == atom2.bondedto or atom == atom2.angleto:
                variablemade[i] = True
                make_variable(atom2)
        global con
        # Next, modify the template connectivity
        # Loop over all connectivities in the template
        for con in tem1.templates[0].connectivity:
            if atom not in con.atoms:
                continue
            atomtypes0 = ["dum" if isinstance(catom.param0, int) else catom.param0.params[
                0] for catom in con.atoms]
            atomtypes1 = ["dum" if isinstance(catom.param1, int) else catom.param1.params[
                0] for catom in con.atoms]
            if isinstance(atom.param1, int):  # Check if we have inserted a dummy parameter
                if con.type == 'bond':
                    con.param0 = con.param1 = 0  # This connectivity should not be sampled
            else:  # Have parameters in cljparams...
                if "dum" in atomtypes0 or "dum" in atomtypes1:
                    continue  # Take care of this if above
                con.param0, equil0 = find_param(atomtypes0, temsets[con.type], gaffsets[
                                                con.type])  # Find parameter index and equilibrium value
                con.param1, equil1 = find_param(
                    atomtypes1, temsets[con.type], gaffsets[con.type])
                if con.param0 == -1 or con.param1 == -1:  # Warn if it could not be found
                    logger.warning("")
                    logger.warning("Warning: could not find parameters for %s or %s at atoms %s" % (
                        "-".join(atomtypes0), "-".join(atomtypes1), " ".join(catom.name for catom in con.atoms)))
                    con.param0 = con.param1 = None

    return tem1


def _make_comb_tem(eletem, vdwtem):
    """
    Make new single topology template file for combined perturbation

    Parameters
    ----------
    eletem : TemplateFile
      the electrostatic template
    vdwtem : TemplateFile
      the van der Waals template

    Returns
    -------
    TemplateFile
      the created template file
    """

    combtem = copy.deepcopy(vdwtem)
    nlj = len(eletem.cljparams)
    nlj2 = nlj / 2

    # Replaces the first half of the CLJ parameters, i.e. the V0 state
    for i in range(nlj2):
        combtem.cljparams[i].params = copy.deepcopy(eletem.cljparams[i].params)

    # Replaces the V0 param of all the atoms
    for atom_ele, atom_comb in zip(eletem.templates[0].atoms, combtem.templates[0].atoms):
        if isinstance(atom_ele.param0, int):
            atom_comb.param0 = atom_ele.param0
        else:
            for clj in combtem.cljparams:
                if clj.index == atom_ele.param0.index:
                    atom_comb.param0 = clj
                    break

    return combtem


def make_single(tem1, tem2, pdb1, pdb2, mapfile=None, gaffversion="gaff16"):
    """
    Make single topology templates

    Parameters
    ----------
    tem1 : TemplateFile or string
      first template file or its filename
    tem2 : TemplateFile or string
      second template file or its filename
    pdb1 : PDBFile or string
      structure template for tem1 or its filename
    pdb2 : PDBFile or string
      structure template for tem2 or its filename
    mapfile : string, optional
      filename of map of correspondance
    gaffversion : str, optional
      version of gaff to use, defaults to gaff16

    Returns
    -------
    TemplateFile
      a template file for the electrostatic leg
    TemplateFile
      a template file for the van der Waals leg
    TemplateFile
      a template file for a combined transformation
    dictionary of strings
      the full map of correspondance

    Raises
    ------
    SetupError
      tem2 is larger than tem1
    """

    logger.debug("Running make_single with arguments: ")
    logger.debug("\ttem1    = %s" % tem1)
    logger.debug("\ttem2    = %s" % tem2)
    logger.debug("\tpdb1    = %s" % pdb1)
    logger.debug("\tpdb2    = %s" % pdb2)
    logger.debug("\tmapfile = %s" % mapfile)
    logger.debug(
        "This will make a ProtoMS template files for single-topology perturbations")

    if isinstance(tem1, basestring):
        tem1 = sim.TemplateFile(filename=tem1)
    if isinstance(tem2, basestring):
        tem2 = sim.TemplateFile(filename=tem2)

    if len(tem1.templates[0].atoms) < len(tem2.templates[0].atoms):
        msg = "The first template needs to be larger than the second template"
        logger.error(msg)
        raise sim.SetupError(msg)

    if isinstance(pdb1, basestring):
        pdb1 = sim.PDBFile(filename=pdb1)
    if isinstance(pdb2, basestring):
        pdb2 = sim.PDBFile(filename=pdb2)

    cmap = {}
    # Read the map from disc
    if mapfile is not None:
        ndum = 0
        with open(mapfile, "r") as f:
            line = f.readline()
            while line:
                atm1, atm2 = line.strip().split()
                atm1 = atm1.upper()
                atm2 = atm2.upper()
                if atm1.startswith("DUM"):
                    ndum = ndum + 1
                    cmap["DUM%s" % ndum] = atm2
                else:
                    cmap[atm1] = atm2
                line = f.readline()

    _make_map(tem1, tem2, pdb1, pdb2, cmap)
    # eletem = _make_ele_tem(tem1, tem2, cmap)
    # vdwtem = _make_vdw_tem(tem1, tem2, pdb1, pdb2, cmap,
    #                        gaffversion=gaffversion)
    # combtem = _make_comb_tem(eletem, vdwtem)
    templates = {'ele': _make_ele_tem(tem1, tem2, cmap),
                 'vdw': _make_vdw_tem(tem1, tem2, pdb1, pdb2, cmap,
                                      gaffversion=gaffversion)}
    templates['comb'] = _make_comb_tem(templates['ele'], templates['vdw'])
    return templates, cmap


def write_map(cmap, filename):
    """
    Write a correspondance map to disc

    Parameters
    ----------
    cmap : a dictionary of strings
      the map of correspondance
    filename : string
      the name of the file
    """

    with open(filename, "w") as f:
        for atm1 in sorted(cmap.keys()):
            atm2 = cmap[atm1]
            if atm1.startswith("DUM"):
                f.write("DUM %s\n" % atm2)
            else:
                f.write("%s %s\n" % (atm1, atm2))


def summarize(eletem, vdwtem, loggfunc):
    """
    Print single topology parameter summary

    Parameters
    ----------
    eletem : TemplateFile 
      template for electrostatic perturbation
    vdwtem : TemplateFile
      template for van der Waals perturbation
    loggfunc : logger function
      the logger function to use to print the summary
    """
    loggfunc("")
    loggfunc("Atom    Ele perturbation                                        ||         Vdw perturbation")
    for atom1, atom2 in zip(eletem.templates[0].atoms, vdwtem.templates[0].atoms):
        ele0, ele1 = atom1.param0, atom1.param1
        vdw0, vdw1 = atom2.param0, atom2.param1

        def get_param(param):
            if isinstance(param, int):
                return "%7.3f %7.3f %7.3f" % (0.000, 0.000, 0.000)
            else:
                return "%7.3f %7.3f %7.3f" % (float(param.params[2]), float(param.params[3]), float(param.params[4]))

        def get_diff(param1, param2):
            if isinstance(param2, int) and not isinstance(param1, int):
                return "***"
            else:
                strout = ""
                if abs(float(param1.params[2]) - float(param2.params[2])) > 0.000001:
                    strout += "c"
                if abs(float(param1.params[3]) - float(param2.params[3])) > 0.000001:
                    strout += "lj"
                return "%3s" % strout
        loggfunc("%5s : %s ==> %s %s || %s ==> %s %s" % (atom1.name, get_param(ele0), get_param(
            ele1), get_diff(ele0, ele1), get_param(vdw0), get_param(vdw1), get_diff(vdw0, vdw1)))

    loggfunc("")
    loggfunc("Z-matrix connectivities that changes parameter: ")
    for con in vdwtem.templates[0].connectivity:
        if con.param0 is not None:
            loggfunc(con)

    loggfunc("")
    loggfunc("Variable geometries: ")
    for comment, var in zip(vdwtem.templates[0].variables[:-1:2], vdwtem.templates[0].variables[1::2]):
        loggfunc("%s %s" % (var, comment))


def get_arg_parser():
    import argparse
    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to setup template files for single-toplogy perturbations semi-automatically")
    parser.add_argument('-t0', '--tem0', help="Template file for V0")
    parser.add_argument('-t1', '--tem1', help="Template file for V1")
    parser.add_argument('-p0', '--pdb0', help="PDB-file for V0")
    parser.add_argument('-p1', '--pdb1', help="PDB-file for V1")
    parser.add_argument(
        '-m', '--map', help="the correspondance map from V0 to V1")
    parser.add_argument(
        '-o', '--out', help="prefix of the output file", default="single")
    parser.add_argument(
        '--gaff', help="the version of GAFF to use for ligand", default="gaff16")
    return parser
#
# If this is run from the command-line
#
if __name__ == '__main__':

    import copy

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = sim.setup_logger("make_single_py.log")

    if args.tem0 is None or args.tem1 is None or args.pdb0 is None or args.pdb1 is None:
        sim.SetupError(
            "Not all four necessary input files given. Cannot setup single-topology")

    templates, cmap = make_single(
        args.tem0, args.tem1, args.pdb0, args.pdb1, args.map, gaffversion=args.gaff)
    for key in templates:
        templates[key].write("%s_%s.tem" % (args.out, key))
    if args.map is None:
        write_map(cmap, args.out + "_cmap.dat")

    # Write out a summary
    summarize(templates['ele'], templates['vdw'], logger.info)
