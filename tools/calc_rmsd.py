# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
""" 
Routine to calculate the RMSD of the center of a ligand

This module defines a these public functions:
calc_rmsd

Can be executed from the command line as a stand-alone program
"""

import logging
import numpy as np
from protomslib import simulationobjects

logger = logging.getLogger('protoms')


def calc_rmsd(ref, structures, resname, atomname=None):
    """
    Calculate RMSD of ligand centre or a specific atom

    Parameters
    ----------
    ref : PDBFile object
      the reference structure
    structures : PDBSet object
      the trajectory to analyze
    resname : string
      the residue name of the ligand to analyze
    atomname : string, optional
      the atom name to analyze, if not set will analyze geometric centre

    Returns
    -------
    float :
      the RMSD in Angstroms
  """
    resname = resname.lower()
    if atomname is not None:
        atomname = atomname.lower()

    centers = []
    for pdb in structures.pdbs:
        for i, res in pdb.residues.iteritems():
            if res.name.lower() != resname:
                continue
            center = np.zeros(3)
            for atom in res.atoms:
                if atomname is not None and \
                   atom.name.strip().lower() != atomname:
                    continue
                center = center + atom.coords
            if atomname is None:
                center = center / float(len(res.atoms))
            centers.append(center)

    refcent = np.zeros(3)
    for i, res in ref.residues.iteritems():
        if res.name.lower() == resname:
            for atom in res.atoms:
                if atomname is not None and \
                   atom.name.strip().lower() != atomname:
                    continue
                refcent = refcent + atom.coords
            if atomname is None:
                refcent = refcent / float(len(res.atoms))
            break

    diff2 = np.sum((centers - refcent)**2, axis=1)
    return np.sqrt(diff2.mean())


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to calculate RMSD of ligand centre")
    parser.add_argument(
        '-i', '--initial', help="the initial PDB-file of the ligand")
    parser.add_argument('-f', '--files', nargs="+", help="the input PDB-files")
    parser.add_argument(
        '-l', '--ligand', help="the name of the ligand to extract")
    parser.add_argument('-a', '--atom', help="the name of the atom to analyze")
    parser.add_argument(
        '-t',
        '--temperature',
        type=float,
        help="the temperature in the simulation",
        default=298.0)
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    if not args.files or not args.initial:
        print("No input files! Nothing to do, so exit.")
        quit()

    if not args.ligand:
        print("No residue name specified, so exit.")
        quit()

    # Read in PDB files
    if len(args.files) == 1:
        pdbfiles = simulationobjects.PDBSet()
        pdbfiles.read(args.files[0], resname=args.ligand)
    else:
        pdbfiles = simulationobjects.PDBSet()
        for filename in args.files:
            pdb = simulationobjects.PDBFile(filename=filename)
            pdbfiles.pdbs.append(pdb)

    # Read in initial structure
    initpdb = simulationobjects.PDBFile(filename=args.initial)

    # Calculate RMSD
    rmsd = calc_rmsd(initpdb, pdbfiles, args.ligand, args.atom)
    if args.atom is None:
        print("The RMSD of the ligand centre is %0.3f A" % rmsd)
    else:
        print("The RMSD of the atom %s is %0.3f A" % (args.atom.lower(), rmsd))

    # Calculate force constant from equipartition theorem
    k = 3.0 * 1.987 * args.temperature / 1000.0 / (rmsd * rmsd)
    print("\nThis corresponds to a spring constant of %.3f kcal/mol/A2" % k)
    print("(to use this as an harmonic restraint in ProtoMS, specify %.3f)" %
          (k * 0.5))
