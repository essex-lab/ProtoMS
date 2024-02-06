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

from protomslib import simulationobjects
from protomslib.prepare import make_dummy

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program make a dummy corresponding to a molecule"
    )
    parser.add_argument("-f", "--file", help="the name of a PDB file")
    parser.add_argument(
        "-o",
        "--out",
        help="the name of the dummy PDB file",
        default="dummy.pdb",
    )
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("make_dummy_py.log")

    if args.file is None:
        print("Nothing to do! Exiting.")

    dummy = make_dummy(args.file)
    dummy.write(args.out)
