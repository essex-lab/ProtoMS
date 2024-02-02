#!/usr/bin/env python
# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to make a GCMC/JAWS-1 simulation box from a ligand

This module defines a single public function:
make_gcmcbox

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np

from protomslib import simulationobjects
from protomslib.gcmc import make_gcmcbox, print_bequil

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to make a PDB-file with box coordinates "
        "covering a solute molecules"
    )
    parser.add_argument(
        "-s",
        "--solute",
        help="the name of the PDB-file containing the solute.",
    )
    parser.add_argument(
        "-p",
        "--padding",
        type=float,
        help="the padding in A,default=2",
        default=2.0,
    )
    parser.add_argument(
        "-o",
        "--out",
        help="the name of the box PDB-file",
        default="gcmc_box.pdb",
    )
    parser.add_argument(
        "-b",
        "--box",
        nargs="+",
        help="Either the centre of the box (x,y,z), or the centre of box AND "
        "length (x,y,z,x,y,z). If the centre is specified and the length "
        "isn't, twice the 'padding' will be the lengths of a cubic box.",
        default=None,
    )
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("make_gcmcbox_py.log")

    if args.solute is None and args.box is None:
        print("Nothing to do! Exiting.")

    if args.box is None:
        pdbobj = simulationobjects.PDBFile()
        pdbobj.read(args.solute)
        box = make_gcmcbox(pdbobj, args.out, args.padding)
    elif len(args.box) == 3:
        box = {
            "center": np.array(
                [float(args.box[0]), float(args.box[1]), float(args.box[2])]
            ),
            "len": np.array([args.padding * 2] * 3),
        }
    elif len(args.box) == 6:
        box = {
            "center": np.array(
                [float(args.box[0]), float(args.box[1]), float(args.box[2])]
            ),
            "len": np.array(
                [float(args.box[3]), float(args.box[4]), float(args.box[5])]
            ),
        }
    else:
        print(
            "\nError with 'box' arguement. Please specify either three "
            "arguments for the centre of the box, or six arguements for "
            "the centre of the box AND the lengths of the sides.\n"
        )

    # Save it to disc
    simulationobjects.write_box(args.out, box)
    print_bequil(box["len"])
