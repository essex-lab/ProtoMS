# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to remove solvent molecules from a GCMC/JAWS-1 simulation box

This module defines a single public function:
clear_gcmcbox

Can be executed from the command line as a stand-alone program
"""

import logging

from protomslib import simulationobjects
from protomslib.prepare.gcmc import clear_gcmcbox

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to remove water molecules from a GCMC/JAWS-1 box"
    )
    parser.add_argument(
        "-b", "--box", help="the name of the PDB-file containing the box."
    )
    parser.add_argument(
        "-s",
        "--solvation",
        help="the name of the PDB-file containing the solvation waters",
    )
    parser.add_argument(
        "-o",
        "--out",
        help="the name of the output PDB-file",
        default="cleared_box.pdb",
    )
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("clear_gcmcbox_py.log")

    if args.solvation is None:
        print("No pdb with solvent provided. Nothing to do! Exiting.")

    nrem, cleared_box = clear_gcmcbox(args.box, args.solvation)
    cleared_box.write(args.out)
