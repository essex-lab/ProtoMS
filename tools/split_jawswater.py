# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to split JAWS-1 water into PDBs for JAWS-2

This module defines a single public function:
split_water

Can be executed from the command line as a stand-alone program
"""

import logging

from protomslib import simulationobjects
from protomslib.prepare.water import split_waters, set_jaws2_box

logger = logging.getLogger('protoms')


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to split JAWS-1 waters to a number of "
                    "PDB-files for JAWS-2")
    parser.add_argument(
        '-w',
        '--waters',
        help="the name of the PDB-file containing the waters.")
    parser.add_argument(
        '-o', '--out', help="the prefix of the output PDB-files", default="")
    parser.add_argument(
        '--jaws2box',
        action='store_true',
        help="whether to apply a header box for jaws2 to the pdb "
             "files of individual waters",
        default=False)
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("split_jawswater_py.log")

    if args.waters is None:
        print("Nothing to do! Exiting.")

    single_waters, other_waters = split_waters(args.waters)

    if args.jaws2box:
        set_jaws2_box(single_waters)

    single_waters.write([
        args.out + "wat%d.pdb" % (i + 1)
        for i in range(len(single_waters.pdbs))
    ])
    other_waters.write([
        args.out + "not%d.pdb" % (i + 1)
        for i in range(len(single_waters.pdbs))
    ])
