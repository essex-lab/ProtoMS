from __future__ import division, print_function

import copy
import logging
import os

import numpy as np
import six

from . import simulationobjects as sim

logger = logging.getLogger("protoms")


def merge_templates(templates):
    """
    Merge a series of ProtoMS template files

    Parameters
    ----------
    templates : list of string
      the names of the template file

    Returns
    -------
    TemplateFile
      the merged template file
    or
      string
    the filename of the single unique template
    """

    logger.debug("Running merge_templates with arguments: ")
    logger.debug("\ttemplates = %s" % " ".join(templates))
    logger.debug(
        "This will merge all templates, renumbering force field parameters"
    )

    # Make it a unique list
    templates2 = []
    for t in templates:
        if t not in templates2:
            templates2.append(t)

    if len(templates2) == 1:
        return templates2[0]

    temfile = sim.TemplateFile(templates2[0])
    for t in templates2[1:]:
        temfile2 = sim.TemplateFile(t)
        temfile.append(temfile2)
    return temfile


class PrepiAtom:
    """
    Class to encapsulate an atom in an Amber prepi file
    and its connectivity

    Attributes
    ----------
    atype : string
      the atom type
    bondidx : int
      the index of the atom this atom is bonded to in the prepi z-matrix
    bonds : list
      names of all atoms this atom is bonded to
    charge : float
      the charge
    index : int
      the serial number in the prepi file
    is_on_loop : dictionary
      a flag for each atom in bonds indicating if the bond is part of a loop
    loop_closure : dictionary
      a flag for each atom in bonds indicating if the bond was part of a loop
      statement in the prepi file
    name : string
      the atom name
    traversed : dictionary
      a flag for each atom in bonds indicating if the bond has been traversed
    zmat : list
      atom names of the atoms defined to this atom in the z-matrix
    """

    def __init__(self):
        self.charge = 0.0
        self.name = ""
        self.index = 0
        self.atype = ""
        self.bondidx = -1
        self.bonds = []
        self.loop_closure = {}
        self.traversed = {}
        self.zmat = []
        self.is_on_loop = {}

    def read(self, record, atomlist):
        """Read a line from an Amber prepi-file"""
        cols = record.strip().split()
        self.name = cols[1]
        self.index = int(cols[0])
        self.atype = cols[2]
        self.bondidx = int(cols[4])
        if self.bondidx >= 4:
            self.add_bond(atomlist[self.bondidx - 1])
        self.charge = float(cols[10])

    def add_bond(self, atomname, loop_closure=False):
        """Add a bond to this atom"""
        self.bonds.append(atomname)
        self.loop_closure[atomname] = loop_closure
        self.is_on_loop[atomname] = False
        self.traversed[atomname] = False

    def sort_bonds(self, metric):
        """Sort the bonds based on some metric"""
        self.bonds = sorted(self.bonds, key=lambda x: metric[x])[::-1]

    def next_bond(self, defined, update=True):
        """Find an atom, bonded to this atom that
        has not been defined and where the bond has not been traversed
        """
        next = None
        for bond in self.bonds:
            if (
                not defined[bond]
                and not self.traversed[bond]
                and self.is_on_loop[bond]
            ):
                if update:
                    self.traversed[bond] = True
                next = bond
                break

        if next is None:
            for bond in self.bonds:
                if (
                    not defined[bond]
                    and not self.traversed[bond]
                    and not self.is_on_loop[bond]
                ):
                    if update:
                        self.traversed[bond] = True
                    next = bond
                    break
        return next

    def backward_bond(self, exclude, defined):
        """
        Find an atom, bonded to this atom that has
        already been defined but is not in the list of excluded atoms
        """
        back = None
        for bond in self.bonds:
            if bond not in exclude and defined[bond]:
                back = bond
                break
        return back

    def improper_dihedral(self, exclude, defined):
        """
        Look for the definition of an improper dihedral:
        see if at least two bonds to this atom has been defined,
        excluding the bond in the exclude list
        """
        a1 = None
        a2 = None
        for bond in self.bonds:
            if bond in exclude:
                continue
            if defined[bond]:
                if a1 is None:
                    a1 = bond
                elif a2 is None:
                    a2 = bond
                    return (a1, a2)
        return (a1, a2)


def _find_cycles(atoms):

    # Create a simple, unweighted, molecular graph
    graph = {}
    for atom in atoms:
        graph[atom] = {}
        for bond in atoms[atom].bonds:
            graph[atom][bond] = {}

    gnodes = set(graph)
    cycles = []
    root = None
    while gnodes:  # loop over connected components
        if root is None:
            root = gnodes.pop()
        stack = [root]
        pred = {root: root}
        used = {root: set()}
        while stack:  # walk the spanning tree finding cycles
            z = stack.pop()  # use last-in so cycles easier to find
            zused = used[z]
            for nbr in graph[z]:
                if nbr not in used:  # new node
                    pred[nbr] = z
                    stack.append(nbr)
                    used[nbr] = set([z])
                elif nbr == z:  # self loops
                    cycles.append([z])
                elif nbr not in zused:  # found a cycle
                    pn = used[nbr]
                    cycle = [nbr, z]
                    p = pred[z]
                    while p not in pn:
                        cycle.append(p)
                        p = pred[p]
                    cycle.append(p)
                    cycles.append(cycle)
                    used[nbr].add(z)
        gnodes -= set(pred)
        root = None
    return cycles


def _read_prepi(filename):
    """
    Read an prepi-file from Antechamber
    """

    atomlist = "DUM DUM DUM".split()
    atoms = {}
    impropers = []

    lines = open(filename, "r").readlines()
    i = 0
    while lines[i].find("CORRECT     OMIT") == -1:
        i = i + 1
    i = i + 5

    while len(lines[i]) > 4:
        atom = PrepiAtom()
        atom.read(lines[i], atomlist)
        atoms[atom.name] = atom
        atomlist.append(atom.name)
        if atom.bondidx >= 4:
            atoms[atom.bonds[-1]].add_bond(atom.name)
        i = i + 1

    while lines[i].find("LOOP") == -1:
        i = i + 1
    i = i + 1
    while len(lines[i]) > 4:
        atom1, atom2 = lines[i].strip().split()
        atoms[atom1].add_bond(atom2, loop_closure=True)
        atoms[atom2].add_bond(atom1)
        i = i + 1

    # Find cycles in the atom graph
    cycles = _find_cycles(atoms)
    for cycle in cycles:
        for atom1, atom2 in zip(cycle[:-1], cycle[1:]):
            atoms[atom1].is_on_loop[atom2] = True
            atoms[atom2].is_on_loop[atom1] = True
        atoms[cycle[0]].is_on_loop[cycle[-1]] = True
        atoms[cycle[-1]].is_on_loop[cycle[0]] = True

    while lines[i].find("IMPROPER") == -1:
        i += 1
    while True:
        i += 1
        cols = lines[i].split()
        if len(cols) != 4:
            break

        impropers.append(cols)

    return atomlist, atoms, impropers


#
# Calculate closeness centrality of all atoms in the molecular graph
#
def _closeness_centrality(graph):
    # Return the length of the shortest path from node to all other nodes
    def shortest_path(graph, node):
        seen = {}
        level = 0
        nextlevel = {node: 1}
        while nextlevel:
            thislevel = nextlevel
            nextlevel = {}
            for v in thislevel:
                if v not in seen:
                    seen[v] = level
                    nextlevel.update(graph[v])
            level = level + 1
        return seen

    c = {}
    for node in graph:
        sp = shortest_path(graph, node)
        spsum = sum(sp.values())
        c[node] = (len(sp) - 1.0) / spsum

    return c


#
# Compute the closeness of all atoms, sort
# their bonds on closeness and return of sorted list of atoms
#
def _compute_closeness(atoms, verbose=False):

    # Create a simple, unweighted, molecular graph
    graph = {}
    for atom in atoms:
        graph[atom] = {}
        for bond in atoms[atom].bonds:
            graph[atom][bond] = {}

    # Calculate the closeness of all atoms
    closeness = _closeness_centrality(graph)

    # Create a sorted list based on the closeness of all atoms
    # sort on atom names as secondary keys to make the resulting
    # order more deterministic, useful for testing
    atomlist = sorted(closeness, key=lambda x: (closeness[x], x))[::-1]

    # Sort the bonds for each atom by the closeness and set
    # all atoms to undefined
    for atom in atomlist:
        atoms[atom].sort_bonds(closeness)
        if verbose:
            print("%s %.3f" % (atom, closeness[atom]))
    if verbose:
        print()

    return atomlist


def make_zmat(prepifile):
    """Make a ProtoMS z-matrix for a solute

    Parameters
    ----------
    prepifile : string
      the filename of the Amber prepi file (prep-file with z-matrix)

    Returns
    -------
    atoms : dictionary
      created PrepiAtom objects
    all_zmat : list
      a list of atoms with their z-matrix atoms
    """

    # ----------------
    # Helper routines
    # ----------------

    def define_atom(
        current, previous, atoms, defined, all_zmat, terminal, verbose=False
    ):
        """Define a new atom in the z-matrix"""

        # Make a proper dihedral for the current atom
        def make_proper():
            # If we a dealing with the first three atoms, just
            # include less and less dummies
            if len(all_zmat) < 2 or (len(all_zmat) == 2 and not terminal):
                atoms[current].zmat = atoms[previous].zmat[:2]
                atoms[current].zmat.insert(0, previous)
            # if not, we traverse backwards on the molecular graph, looking
            # for two atoms to form an angle and dihedral
            else:
                a2 = atoms[previous].backward_bond(
                    [current, previous], defined
                )
                if len(all_zmat) > 2:
                    a3 = atoms[a2].backward_bond([previous, a2], defined)
                else:
                    a3 = "DM3"
                atoms[current].zmat = ("%s %s %s" % (previous, a2, a3)).split()

        # Make an improper dihedral for the current atom
        def make_improper(a2, a3):
            # This two lines were the previous default definition,
            # now they are supplied as arguments
            atoms[current].zmat = ("%s %s %s" % (previous, a2, a3)).split()

        defined[current] = True
        if previous is not None:

            if len(all_zmat) <= 2:
                make_proper()
            else:

                # This was the previous rule:
                # If the current atom is the second or more bonded atom to
                # the previous, let's define it with an improper so that
                # it moves with the first atom bonded to the previous atom
                # if atoms[previous].bonds.index(current) > 1 :

                # Now, we let the previous atom make the decision:
                # simply return two previously defined atoms, or None
                a2, a3 = atoms[previous].improper_dihedral(current, defined)
                if a2 is not None and a3 is not None:
                    make_improper(a2, a3)
                else:
                    make_proper()

        # This will only happens for the very first atom
        else:
            atoms[current].zmat = "DM3 DM2 DM1".split()
        if verbose:
            print(" ", current, end="")
        all_zmat.append(current + " " + " ".join(atoms[current].zmat))

    def traverse_graph(atomlist, atoms, verbose=False):
        """
        Traverse the molecule graph based on the closeness of the atoms
        and define the z-matrix
        """
        defined = {}
        for atom in atoms:
            defined[atom] = False
        all_zmat = []

        if verbose:
            print("Traversal of the molecular graph:")

        # Start with the most central atom
        branch_atom = atomlist[0]
        terminal_flag = False
        define_atom(
            branch_atom, None, atoms, defined, all_zmat, terminal_flag, verbose
        )

        while True:
            # Traverse a branch from an atom with at least two bonds
            previous = branch_atom
            next = atoms[branch_atom].next_bond(defined)
            atoms[next].traversed[branch_atom] = True
            while next is not None:
                define_atom(
                    next,
                    previous,
                    atoms,
                    defined,
                    all_zmat,
                    terminal_flag,
                    verbose,
                )
                previous = next
                next = atoms[previous].next_bond(defined)
                if next is not None:
                    atoms[next].traversed[previous] = True
                terminal_flag = next is None
            if verbose:
                print(".", end="")

            # Check if we have more branches to traverse
            if atoms[branch_atom].next_bond(defined, update=False) is None:
                if verbose:
                    print(" : ", end="")
                # Tries to find a new atom to branch off from
                found = False
                for atom in atomlist:
                    if not defined[atom]:
                        branch_atom = atom
                        terminal_flag = False
                        define_atom(
                            branch_atom,
                            atoms[atom].bonds[0],
                            atoms,
                            defined,
                            all_zmat,
                            terminal_flag,
                            verbose,
                        )
                        # If the found atom has more than one bonds,
                        # we can traverse its branches
                        if len(atoms[branch_atom].bonds) > 1:
                            found = True
                            break
                        else:
                            if verbose:
                                print("; ", end="")
                # If we could not find any more atoms to branch off
                # from we are done!
                if not found:
                    break
        if verbose:
            print()

        return all_zmat

    logger.debug("Running make_zmat with arguments: ")
    logger.debug("\tprepifile = %s" % prepifile)
    logger.debug(
        "This will generate a ProtoMS compatible z-matrix for a solute"
    )

    # Parse the prepi-file into a list of atom names and
    # a dictionary of Atom objects
    atomnames, atoms, impropers = _read_prepi(prepifile)  # MOD
    h, t = os.path.splitext(prepifile)

    # Compute closeness of all atoms and sort their bonds based on this
    atomlist_closeness = _compute_closeness(atoms, verbose=False)  # MOD
    # Traverse the molecular graph and define the z-matrix
    all_zmat = traverse_graph(atomlist_closeness, atoms, verbose=False)  # MOD

    return atoms, all_zmat, impropers


def _readfrcmod(filename):
    """Read an Amber frcmod file from disc

    readfrcmod(filename)

    Parameters
    ----------
    filename - the filename of the Amber frcmod file

    Returns
    -------
    a list of bond, angle and dihedral parameters
    """

    with open(filename) as f:
        # Find start of bond region
        while not next(f).startswith("BOND"):
            pass

        bonds = []
        for line in f:
            cols = line.split()
            if not cols:
                break
            bonds.append([cols[0].split("-")] + list(map(float, cols[1:3])))
        next(f)

        angles = []
        for line in f:
            cols = line.replace(" -", "-").split()
            if not cols:
                break
            angles.append([cols[0].split("-")] + list(map(float, cols[1:3])))

        next(f)

        dihedrals = []
        for line in f:
            cols = line.replace(" -", "-").split()
            if not cols:
                break
            div, k, n, phi = map(float, cols[1:5])
            dihedrals.append([cols[0].split("-"), 1.0, k / div, n, phi])

        next(f)

    return bonds, angles, dihedrals


def build_template(
    temfile,
    prepifile,
    translate=0.1,
    rotate=1.0,
    zmatfile=None,
    frcmodfile=None,
    resname="UNK",
    alldihs=False,
    gaffversion="gaff16",
):
    """Build a ProtoMS template file

    Parameters
    ----------
    temfile : string
      the filename to save the template file to
    prepifile : string
      the filename of the Amber prepi file (prep-file with z-matrix)
    translate : float, optional
      the translational displacement
    rotate : float, optional
      the rotational displacement
    zmatfile : string, optional
      the filename of a zmat, if None it is created
    frcmodfile : string, optional
      the filename of an Amber frcmod file with additional parameters
    resname : string, optional
      the name of solute
    alldihs : boolean, optional
      set True to sample improper dihedrals
    gaffversion : string, optional
      the version of GAFF to use

    Returns
    -------
    TemplateFile
      the created template file
    """

    logger.debug("Running build_template with arguments: ")
    logger.debug("\ttemfile     = %s" % temfile)
    logger.debug("\tprepifile   = %s" % prepifile)
    logger.debug("\ttranslate   = %f" % translate)
    logger.debug("\trotate      = %f" % rotate)
    logger.debug("\tzmatfile    = %s" % zmatfile)
    logger.debug("\tfrcmodfile  = %s" % frcmodfile)
    logger.debug("\tresname     = %s" % resname)
    logger.debug("\tgaffversion = %s" % gaffversion)
    logger.debug("This will generate a ProtoMS template file for a solute")

    if zmatfile is None:
        atoms, zmat, impropers = make_zmat(prepifile)
    else:
        with open(zmatfile) as f:
            zmat = f.readlines()
        atoms, impropers = _read_prepi(prepifile)[1:]

    centralities = _compute_closeness(atoms, verbose=False)

    if frcmodfile is None:
        frcbonds = frcangles = frcdihedrals = None
    else:
        frcbonds, frcangles, frcdihedrals = _readfrcmod(frcmodfile)

    gaff_file = sim.standard_filename(gaffversion + ".ff", "parameter")
    angle_params = sim.ParameterSet("angle", gaff_file)
    dihedral_params = sim.ParameterSet("dihedral", gaff_file)

    if gaffversion == "gaff14":
        gaffversion = "gaff"
    with open(sim.standard_filename(gaffversion + ".types", "parameter")) as f:
        at_params = [line.split() for line in f]

    kBT = 0.0019872041 * 300  # Boltzmann constant at 300 kelvin in kcal/mol
    move_scale = 0.5

    dummynames = ["DM3", "DM2", "DM1"]
    aromatic_types = (
        "ca cp cq ce cf cc cd nb ne nf pb pe pf px py sx sy"
        "CA CB CC CK CM CN CQ CR CV CW C* NA NB NC"
    ).split()
    template = sim.TemplateFile()
    template.templates.append(sim.MolTemplate())
    template.templates[0].atomclass = sim.TemplateSoluteAtom
    moltem = template.templates[0]
    moltem.type = "solute"
    moltem.name = resname
    moltem.translate = translate
    moltem.rotate = rotate

    def get_resname(s):
        if s in dummynames:
            return "DUM"
        else:
            return resname

    if frcbonds is not None:
        for i, bond in enumerate(frcbonds, 4500):
            template.bondparams.append(
                sim.ForceFieldParameter(
                    record="par %4d %.3f %.3f\n" % (i, bond[1], bond[2])
                )
            )
        for i, bond in enumerate(frcbonds, 4500):
            template.bondatoms.append(
                sim.AtomSet(
                    record="atm %s %s %4d\n" % (bond[0][0], bond[0][1], i)
                )
            )

    if frcangles is not None:
        for i, angle in enumerate(frcangles, 4700):
            template.angleparams.append(
                sim.ForceFieldParameter(
                    record="par %4d %.3f %.3f\n" % (i, angle[1], angle[2])
                )
            )
        for i, angle in enumerate(frcangles, 4700):
            template.angleatoms.append(
                sim.AtomSet(
                    record="atm %s %s %s %4d\n"
                    % (angle[0][0], angle[0][1], angle[0][2], i)
                )
            )

    if frcdihedrals is not None:
        fmt = (4600, 0.0, 0.0, 0.0, 0.0)
        template.dihedralterms.append(
            sim.ForceFieldParameter(
                record="term  %4d   %.3f   %.3f   %.3f   %.3f\n" % fmt
            )
        )
        for i, di in enumerate(frcdihedrals, 4601):
            fmt = (i, di[2], di[1], di[4], di[3])
            template.dihedralterms.append(
                sim.ForceFieldParameter(
                    record="term  %4d   %.3f   %.3f   %.3f   %.3f\n" % fmt
                )
            )

        # Loop through unique dihedral terms
        done = []
        for i, di in enumerate(frcdihedrals, 4600):
            if di[0] not in done:
                # Get all terms with matching atom groups from dihedral list
                terms = [
                    (j, term)
                    for j, term in enumerate(frcdihedrals, 4601)
                    if term[0] == di[0]
                ]
                # Start par entry
                record = "par  %4d" % i
                # Loop through all 5 terms of full dihedral
                for j in range(5):
                    # Pull out term with values of k3 == j
                    ref = [term for term in terms if int(term[1][4]) == j]
                    if ref:
                        record += "   %4d" % ref[0][0]
                    else:
                        # If no k3 == j give default zero parameters
                        record += "   %4d" % 4600
                template.dihedralparams.append(
                    sim.ForceFieldParameter(record=record)
                )
            done += [di[0]]

        done = []
        for i, di in enumerate(frcdihedrals, 4600):
            if di[0] not in done:
                fmt = (di[0][0], di[0][1], di[0][2], di[0][3], i)
                template.dihedralatoms.append(
                    sim.AtomSet(record="atm %4s %4s %4s %4s %4d\n" % fmt)
                )
            done += [di[0]]

    # Here we need to iterate over the atom order in the z-matrix,
    # so this is modified
    for i, atom in enumerate(
        [atoms[zatms.split()[0]] for zatms in zmat], 3000
    ):
        params = [j for j in at_params if j[0] == atom.atype][0]
        fmt = (
            i,
            atom.atype,
            int(params[6]),
            atom.charge,
            float(params[2]),
            float(params[4]),
        )
        template.cljparams.append(
            sim.ForceFieldParameter(
                record="par  %4d %4s   %02d %10.5f %10.5f %10.5f\n" % fmt
            )
        )

    logger.info(
        "Before running a simulation, ensure that the first line of %s.pdb "
        "reads 'HEADER %s'." % (os.path.splitext(temfile)[0], resname)
    )

    # Print out the atoms
    for i, line in enumerate(zmat, 3000):
        atms = line.split()
        fmt = (
            atms[0],
            resname,
            i,
            i,
            atms[1],
            get_resname(atms[1]),
            atms[2],
            get_resname(atms[2]),
            atms[3],
            get_resname(atms[3]),
        )
        moltem.atoms.append(
            sim.TemplateSoluteAtom(
                record="atom  %4s %4s %4i %4i %4s %4s %4s %4s %4s %4s\n" % fmt
            )
        )

    # Print out the bonds
    taken_bonds = []
    for line in zmat[1:]:
        atom = line.split()[:2]
        taken_bonds.append((atom[0], atom[1]))
        taken_bonds.append((atom[1], atom[0]))
        moltem.connectivity.append(
            sim.TemplateConnectivity(
                record="bond %4s %3s %4s %3s\n"
                % (atom[0], resname, atom[1], resname)
            )
        )

    # Also need to do bonds that close loops
    for atom1 in sorted(atoms):
        for atom2 in sorted(atoms[atom1].bonds):
            if (atom1, atom2) not in taken_bonds:
                taken_bonds.append((atom1, atom2))
                taken_bonds.append((atom2, atom1))
                # determine atom order by name, this improves consistency
                # of templates between python versions for testing
                at1, at2 = sorted((atom1, atom2))
                moltem.connectivity.append(
                    sim.TemplateConnectivity(
                        record="bond %4s %3s %4s %3s\n"
                        % (at1, resname, at2, resname)
                    )
                )

    # Print out the angles:
    for line in zmat[2:]:
        angle = line.split()[:3]
        if "DM3" in angle:
            continue

        # if this is a loop angle then exclude it
        if (
            atoms[angle[1]].is_on_loop[angle[0]]
            and atoms[angle[1]].is_on_loop[angle[2]]
        ):
            continue
        atypes = [atoms[i].atype for i in angle]
        try:
            k = angle_params.get_params(atypes).k
        except KeyError:
            # parameter not in gaff.ff try frcmod params
            try:
                k = [i for i in frcangles if atypes in (i[0], i[0][::-1])][0][
                    1
                ]
            except IndexError:
                logger.warning(
                    "Unable to find angle parameters for %s-%s-%s"
                    % tuple(atypes)
                )
                logger.warning(
                    "For now the flexibility of this angle will be set to zero"
                )
                logger.warning(
                    "To correct this consider manually adding a term to the "
                    "frcmod file\n"
                )
                k = 10**10  # make k huge so max_move is zero
            except TypeError:
                k = 10 * 10

        max_move = move_scale * 2 * (1 / (kBT * k)) ** 0.5
        if False in [i in aromatic_types for i in atypes]:
            fmt = (
                angle[0],
                resname,
                angle[1],
                resname,
                angle[2],
                resname,
                max_move,
            )
            moltem.connectivity.append(
                sim.TemplateConnectivity(
                    record="angle %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt
                )
            )

    min_flex = 2.0
    max_flex = 10.0
    for line in zmat[3:]:  # First 3 lines contain dummies
        # The first 4 atoms in the z matrix corresponds to the dihedrals.
        dihedral = line.split()[:4]
        atypes = [atoms[i].atype for i in dihedral]
        # Check matching parameters exist in gaff.ff
        missing = False
        try:
            dihedral_params.get_params(atypes)
        except KeyError:
            # If not in gaff.ff then check frcmod
            try:
                [i for i in frcdihedrals if atypes in (i[0], i[0][::-1])]
            except IndexError:
                logger.warning(
                    "Unable to find dihedral parameters for %s-%s-%s-%s"
                    % tuple(atypes)
                )
                logger.warning(
                    "For now the flexibility of this angle will be set to zero"
                )
                logger.warning(
                    "To correct this consider manually adding a term to the"
                    " frcmod file\n"
                )
                k = 10**10  # make k huge so max_move is zero
                missing = True

        skip = False
        # Ensure improper dihedrals are removed, if desired, by checking the
        # last atom in the dihedral is bonded to the previous atom in the list.
        if alldihs:
            # even if we're sampling improper dihedrals we should not sample
            # those that maintain in-plane geometries as ProtoMS does not
            # evaluate these energies. These appear as impropers in the
            # prepi file
            for imp in impropers:
                if (
                    dihedral[0] in imp
                    and dihedral[1] in imp
                    and dihedral[2] in imp
                    and dihedral[3] in imp
                ):
                    skip = True
                    break

        else:
            if dihedral[2] not in atoms[dihedral[3]].bonds:
                skip = True

        if skip:
            continue

        # If central bond is part of a loop
        if atoms[dihedral[1]].is_on_loop[dihedral[2]]:
            # Ensure rotations around aromatic dihedrals are excluded
            if atypes[1] in aromatic_types and atypes[2] in aromatic_types:
                continue

            # If this is a dihedral of all loop atoms then exclude it as well
            if (
                atoms[dihedral[0]].is_on_loop[dihedral[1]]
                and atoms[dihedral[2]].is_on_loop[dihedral[3]]
            ):
                continue

            fmt = (
                dihedral[0],
                resname,
                dihedral[1],
                resname,
                dihedral[2],
                resname,
                dihedral[3],
                resname,
                min_flex,
            )
            moltem.connectivity.append(
                sim.TemplateConnectivity(
                    record="dihedral %4s %4s %4s %4s %4s %4s %4s %4s "
                    "flex %.3f\n" % fmt
                )
            )
        else:
            # Not part of loop so base flexibility on centrality
            av_rank = (
                centralities.index(dihedral[1])
                + centralities.index(dihedral[2])
            ) / float(2 * len(centralities))
            flex = (av_rank) * (max_flex - min_flex) + min_flex
            if missing:
                flex = 0.0
            fmt = (
                dihedral[0],
                resname,
                dihedral[1],
                resname,
                dihedral[2],
                resname,
                dihedral[3],
                resname,
                flex,
            )
            moltem.connectivity.append(
                sim.TemplateConnectivity(
                    record="dihedral %4s %4s %4s %4s %4s %4s %4s %4s"
                    " flex %.3f\n" % fmt
                )
            )

    return template


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
                "Warning: Could not find atom %s in the pdb-file, cannot"
                " look for map for this atom" % atom
            )
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

    moltem1 = tem1.templates[0]
    moltem2 = tem2.templates[0]
    objdict1 = _make_dict(
        [at.name for at in moltem1.atoms], moltem1, pdb1.residues[1]
    )
    objdict2 = _make_dict(
        [at.name for at in moltem2.atoms], moltem2, pdb2.residues[1]
    )

    not_taken1, not_taken2 = _auto_map(
        moltem1, moltem2, objdict1, objdict2, cmap
    )

    # Now we have gotten so-far that we have to ask the user for the rest...

    if not_taken1 or not_taken2:

        # Print out useful information
        logger.info("")
        logger.info(
            "These are the un-matched atoms of template 1: %s"
            % " ".join(not_taken1)
        )
        logger.info(
            "These are the un-matched atoms of template 2: %s"
            % " ".join(not_taken2)
        )

        logger.info("")
        logger.info("These are the distances (A): ")
        logger.info(
            "%8s%s" % ("", "".join("%8s" % atom for atom in not_taken1))
        )
        for atom2 in not_taken2:
            outstr = "%8s" % atom2
            for atom1 in not_taken1:
                dist = (
                    objdict1[atom1]["pdb"].coords
                    - objdict2[atom2]["pdb"].coords
                )
                dist = np.sqrt(np.sum(dist * dist))
                outstr = outstr + "%8.3f" % dist
            logger.info(outstr)
        logger.info("")

        found = {}
        for atom1 in not_taken1:
            atom2 = ""
            dummy = False
            while atom2 not in not_taken2 and not dummy:
                print("Enter the corresponding atom for %s: " % atom1, end="")
                test = six.moves.input().strip().upper()
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
        logger.info(
            "User maps:\n%s"
            % " ".join(
                "%s-%s" % (atom1, found[atom1])
                for atom1 in sorted(found.keys())
            )
        )
        for atom1 in found:
            atom2 = found[atom1]
            not_taken1.remove(atom1)
            if atom2 != "DUM":
                not_taken2.remove(atom2)
            cmap[atom1] = atom2


def _auto_map(tem1, tem2, objdict1, objdict2, cmap):
    """
    As far as possible update cmap based on atom matching
    between templates and pbs according to the below criteria:
    1) squared inter-atom distances < 0.02 A**2
    2) atom type

    The cmap dictionary will be populated with key value pairs:
    cmap[atom1] = atom2,
    where atom1 belongs to tem1 and atom2 belongs to tem2
    atom2 can be dum, in case non-correspondance is assumed

    Parameters
    ----------
    tem1 : MolTemplate
      first template
    tem2 : MolTemplate
      second template file
    pdb1 : PDBFile
      structure template for tem1
    pdb2 : PDBFile
      structure template for tem2
    cmap : dictionary of string
      full, partial or empty map of correspondance

    Returns
    -------
    not_taken1 : list of strings
      list of unpaired atom names from ligand 1
    not_taken2 : list of strings
      list of unpaired atom names from ligand 2

    Raises
    ------
    SetupError
      if dummies in cmap.keys()
    """
    # Create a list of atom names in both templates
    # The purpose of the routine is then to empty these lists!
    not_taken1 = [atom.name.strip().upper() for atom in tem1.atoms]
    not_taken2 = [atom.name.strip().upper() for atom in tem2.atoms]

    # First we will look in the given cmap to see if we already
    # have some maps given
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
                        "Warning: Could not find atom %s in the second "
                        "template file, will ignore this map." % atom2
                    )
                    del cmap[atom1]
            else:
                not_taken1.remove(atom1)
        else:
            logger.warning(
                "Warning: Could not find atom %s in the first template file, "
                "will ignore this map." % atom2
            )
            del cmap[atom1]
    logger.info(
        "Pre-defined maps:\n%s"
        % " ".join(
            "%s-%s" % (atom1, cmap[atom1]) for atom1 in sorted(cmap.keys())
        )
    )

    # Next we will calculate pair-wise distance from the given pdb structures,
    # use a cut-off of 0.02 A^2 to determine possible pair and
    # then use atom type as the final filter
    if not_taken1 or not_taken2:
        found = {}
        for atom1 in not_taken1:
            if objdict1[atom1]["pdb"] is None:
                continue
            for atom2 in not_taken2:
                if objdict2[atom2]["pdb"] is None:
                    continue
                dist2 = (
                    objdict1[atom1]["pdb"].coords
                    - objdict2[atom2]["pdb"].coords
                )
                dist2 = np.sum(dist2 * dist2)
                # The atom type is the 0:th parameter of a clj parameter
                type1 = objdict1[atom1]["tem"].param0.params[0]
                type2 = objdict2[atom2]["tem"].param0.params[0]
                if dist2 < 0.02 and type1 == type2:
                    found[atom1] = atom2
                    break
        logger.info(
            "Distance/atom type maps:\n%s"
            % " ".join(
                "%s-%s" % (atom1, found[atom1])
                for atom1 in sorted(found.keys())
            )
        )
        for atom1 in found:
            atom2 = found[atom1]
            not_taken1.remove(atom1)
            not_taken2.remove(atom2)
            cmap[atom1] = atom2
    return not_taken1, not_taken2


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

    # Make copies of the objects first so we can manipulate them
    # without modifying the original templates
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
                    atom1.param1.params = list(
                        atom1.param0.params
                    )  # This make sure that we do not change vdw
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


def _make_vdw_tem(
    tem1, tem2, pdb1, pdb2, cmap, usepdb=True, gaffversion="gaff16"
):
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
        """Find force field index and equilibrium value for atom types"""
        ratoms = atoms[::-1]
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
        """Find bond lengths and angles in pdb files"""
        allcoords = [objdict[atom]["pdb"].coords for atom in atoms]
        allcoords = np.array(allcoords)
        if len(atoms) == 2:
            return np.sqrt(np.sum((allcoords[0, :] - allcoords[1, :]) ** 2))
        elif len(atoms) == 3:
            return (
                sim.angle_atms(
                    allcoords[0, :], allcoords[1, :], allcoords[2, :]
                )
                * 180.0
                / np.pi
            )
        return 0.0

    def make_variable(atom):
        """Add variable geometry to a molecular template"""
        # Skip this if the atom is bonded to a dummy atom
        if atom.bondedto in ["DM1", "DM2", "DM3"]:
            return
        # This flags the creating of a dummy angle
        defangle = atom.angleto not in ["DM1", "DM2", "DM3"]

        # Store away z-matrix atoms and their names
        zmatatoms = [atom, atom.bondedto]
        if defangle:
            zmatatoms.append(atom.angleto)
        zmatnames1 = [zatom.name.upper() for zatom in zmatatoms]
        zmatnames2 = [cmap[zatom.name.upper()] for zatom in zmatatoms]
        atomtypes1 = [
            zatom.param0.params[0] for zatom in zmatatoms
        ]  # Atom types for template 1
        if usepdb:  # If we should take the geometry from the pdb-files
            _make_dict(
                zmatnames1,
                tem1.templates[0],
                pdb1.residues[1],
                objdict1,
                onlypdb=True,
            )  # Find the Atom object and put in a dictionary
            bond1 = find_pdbparam(zmatnames1[:2], objdict1)
        else:  # If we should take the geometry from GAFF
            dummy, bond1 = find_param(
                atomtypes1[:2], temsets["bond"], gaffsets["bond"]
            )
        if isinstance(
            atom.param1, int
        ):  # Check if this atom should be perturbed to a dummy
            tem1.templates[0].variables.append(
                "# %s to dummy" % "-".join(atomtypes0[:2])
            )
            tem1.templates[0].variables.append(
                "variable %s %s bond %.3f %.3f"
                % (atom.name, atom.residue, bond1, 0.200)
            )  # Shrink bond to within vdw-sphere
        else:  # State V1 has cljparameters
            atomtypes2 = [zatom.param1.params[0] for zatom in zmatatoms]
            if usepdb:
                _make_dict(
                    zmatnames2,
                    tem2.templates[0],
                    pdb2.residues[1],
                    objdict2,
                    onlypdb=True,
                )
                bond2 = find_pdbparam(zmatnames2[:2], objdict2)
                if defangle:
                    angle1 = find_pdbparam(zmatnames1, objdict1)
                    angle2 = find_pdbparam(zmatnames2, objdict2)
            else:
                dummy, bond2 = find_param(
                    atomtypes2[:2], temsets["bond"], gaffsets["bond"]
                )
                if defangle:
                    dummy, angle1 = find_param(
                        atomtypes1, temsets["angle"], gaffsets["angle"]
                    )
                    dummy, angle2 = find_param(
                        atomtypes2, temsets["angle"], gaffsets["angle"]
                    )
            # Add variable geometry for bond and angle
            if bond1 != bond2:
                tem1.templates[0].variables.append(
                    "# %s to %s at atoms %s"
                    % (
                        "-".join(atomtypes1[:2]),
                        "-".join(atomtypes2[:2]),
                        "-".join(zmatnames1[:2]),
                    )
                )
                tem1.templates[0].variables.append(
                    "variable %s %s bond %.3f %.3f"
                    % (atom.name, atom.residue, bond1, bond2)
                )
            if defangle and angle1 != angle2:
                tem1.templates[0].variables.append(
                    "# %s to %s at atoms %s"
                    % (
                        "-".join(atomtypes1),
                        "-".join(atomtypes2),
                        "-".join(zmatnames1),
                    )
                )
                tem1.templates[0].variables.append(
                    "variable %s %s angle %.3f %.3f"
                    % (atom.name, atom.residue, angle1, angle2)
                )

    # MAIN routine

    # Make copies of the objects first so we can manipulate them
    # without modifying the original templates
    tem1 = copy.deepcopy(tem1)
    tem2 = copy.deepcopy(tem2)

    # Copy all clj parameters from tem2 to tem1, both of the sets
    # will be modified below
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
                    #  Change the charges
                    atom1.param0.params[2] = atom3.param0.params[2]
                    # Change the param1, charge AND vdw
                    atom1.param1 = atom3.param0
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

    # Now we need to change the connectivity of the template so that the
    # intramolecular energies/geometry are correctly calculated

    # Setup parameter sets to search
    gaffname = sim.standard_filename(gaffversion + ".ff", "parameter")
    gaffsets = {
        ptype: sim.ParameterSet(ptype, gaffname)
        for ptype in ["bond", "angle", "dihedral"]
    }
    temsets = {
        "bond": tem1.bondatoms,
        "angle": tem1.angleatoms,
        "dihedral": tem1.dihedralatoms,
    }

    # Now loop over all atoms that changes type
    objdict1 = _make_dict(
        [atom.name for atom in atomlist],
        tem1.templates[0],
        pdb1.residues[1],
        onlypdb=True,
    )
    objdict2 = _make_dict(
        [
            cmap[atom.name.upper()]
            for atom in atomlist
            if cmap[atom.name.upper()] != "DUM"
        ],
        tem2.templates[0],
        pdb2.residues[1],
        onlypdb=True,
    )
    variablemade = [False] * len(tem1.templates[0].atoms)
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
        for con in tem1.templates[0].connectivity:
            if atom not in con.atoms:
                continue
            atomtypes0 = [
                (
                    "dum"
                    if isinstance(catom.param0, int)
                    else catom.param0.params[0]
                )
                for catom in con.atoms
            ]
            atomtypes1 = [
                (
                    "dum"
                    if isinstance(catom.param1, int)
                    else catom.param1.params[0]
                )
                for catom in con.atoms
            ]
            if isinstance(
                atom.param1, int
            ):  # Check if we have inserted a dummy parameter
                if con.type == "bond":
                    # This connectivity should not be sampled
                    con.param0 = con.param1 = 0
            else:  # Have parameters in cljparams...
                if "dum" in atomtypes0 or "dum" in atomtypes1:
                    continue  # Take care of this if above
                con.param0, equil0 = find_param(
                    atomtypes0, temsets[con.type], gaffsets[con.type]
                )  # Find parameter index and equilibrium value
                con.param1, equil1 = find_param(
                    atomtypes1, temsets[con.type], gaffsets[con.type]
                )
                if con.param0 == -1 or con.param1 == -1:
                    # Warn if it could not be found
                    logger.warning("")
                    logger.warning(
                        "Warning: could not find parameters for %s or %s at"
                        "atoms %s"
                        % (
                            "-".join(atomtypes0),
                            "-".join(atomtypes1),
                            " ".join(catom.name for catom in con.atoms),
                        )
                    )
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
    nlj2 = nlj // 2

    # Replaces the first half of the CLJ parameters, i.e. the V0 state
    for i in range(nlj2):
        combtem.cljparams[i].params = copy.deepcopy(eletem.cljparams[i].params)

    # Replaces the V0 param of all the atoms
    for atom_ele, atom_comb in zip(
        eletem.templates[0].atoms, combtem.templates[0].atoms
    ):
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
        "This will make a ProtoMS template files for "
        "single-topology perturbations"
    )

    if isinstance(tem1, six.string_types):
        tem1 = sim.TemplateFile(filename=tem1)
    if isinstance(tem2, six.string_types):
        tem2 = sim.TemplateFile(filename=tem2)

    if len(tem1.templates[0].atoms) < len(tem2.templates[0].atoms):
        msg = "The first template needs to be larger than the second template"
        logger.error(msg)
        raise sim.SetupError(msg)

    if isinstance(pdb1, six.string_types):
        pdb1 = sim.PDBFile(filename=pdb1)
    if isinstance(pdb2, six.string_types):
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
    eletem = _make_ele_tem(tem1, tem2, cmap)
    vdwtem = _make_vdw_tem(
        tem1, tem2, pdb1, pdb2, cmap, gaffversion=gaffversion
    )
    combtem = _make_comb_tem(eletem, vdwtem)

    return eletem, vdwtem, combtem, cmap


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


def summarize_single(eletem, vdwtem, loggfunc):
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
    loggfunc(
        "Atom    Ele perturbation                                        ||         Vdw perturbation"
    )
    for atom1, atom2 in zip(
        eletem.templates[0].atoms, vdwtem.templates[0].atoms
    ):
        ele0, ele1 = atom1.param0, atom1.param1
        vdw0, vdw1 = atom2.param0, atom2.param1

        def get_param(param):
            if isinstance(param, int):
                return "%7.3f %7.3f %7.3f" % (0.000, 0.000, 0.000)
            else:
                return "%7.3f %7.3f %7.3f" % (
                    float(param.params[2]),
                    float(param.params[3]),
                    float(param.params[4]),
                )

        def get_diff(param1, param2):
            if isinstance(param2, int) and not isinstance(param1, int):
                return "***"
            else:
                strout = ""
                if (
                    abs(float(param1.params[2]) - float(param2.params[2]))
                    > 0.000001
                ):
                    strout += "c"
                if (
                    abs(float(param1.params[3]) - float(param2.params[3]))
                    > 0.000001
                ):
                    strout += "lj"
                return "%3s" % strout

        loggfunc(
            "%5s : %s ==> %s %s || %s ==> %s %s"
            % (
                atom1.name,
                get_param(ele0),
                get_param(ele1),
                get_diff(ele0, ele1),
                get_param(vdw0),
                get_param(vdw1),
                get_diff(vdw0, vdw1),
            )
        )

    loggfunc("")
    loggfunc("Z-matrix connectivities that changes parameter: ")
    for con in vdwtem.templates[0].connectivity:
        if con.param0 is not None:
            loggfunc(con)

    loggfunc("")
    loggfunc("Variable geometries: ")
    for comment, var in zip(
        vdwtem.templates[0].variables[:-1:2],
        vdwtem.templates[0].variables[1::2],
    ):
        loggfunc("%s %s" % (var, comment))
