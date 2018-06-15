# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to make a dummy PDB structure for a solute

This module defines a single public function:
make_dummy

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np

from protomslib import simulationobjects
import six

logger = logging.getLogger('protoms')


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


def get_arg_parser():
    import argparse
    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program make a dummy corresponding to a molecule")
    parser.add_argument('-f', '--file', help="the name of a PDB file")
    parser.add_argument(
        '-o',
        '--out',
        help="the name of the dummy PDB file",
        default="dummy.pdb")
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("make_dummy_py.log")

    if args.file is None:
        print("Nothing to do! Exiting.")

    dummy = make_dummy(args.file)
    dummy.write(args.out)
