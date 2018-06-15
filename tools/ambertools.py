# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines that wrap AmberTools programs

This module defines two public functions:
run_antechamber
run_parmchk

Can be executed from the command line as a stand-alone program
"""

import logging

from protomslib import simulationobjects
from protomslib.prepare import run_antechamber, run_parmchk

logger = logging.getLogger('protoms')


def get_arg_parser():

    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to run antechamber and parmchk"
                    " for a series of PDB-files")
    parser.add_argument(
        '-f', '--files', nargs="+", help="the name of the PDB-files")
    parser.add_argument(
        '-n', '--name', help="the name of the solute", default="UNK")
    parser.add_argument(
        '-c',
        '--charge',
        nargs="+",
        type=float,
        help="the net charge of each PDB-file")
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    # Setup the logger
    logger = simulationobjects.setup_logger("ambertools_py.log")

    for i, filename in enumerate(args.files):
        if args.charge is None or i >= len(args.charge):
            charge = 0.0
        else:
            charge = args.charge[i]
        run_antechamber(filename, charge, args.name)
        run_parmchk(filename + "P")
