# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines to build a protein scoop

This module defines a single public function
scoop

Can be executed from the command line as a stand-alone program
"""
import logging

from protomslib import simulationobjects
from protomslib.prepare import scoop


logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program scoop a protein pdb-file"
    )
    parser.add_argument("-p", "--protein", help="the protein PDB-file")
    parser.add_argument("-l", "--ligand", help="the ligand PDB-file")
    parser.add_argument(
        "-o", "--out", help="the output PDB-file", default="scoop.pdb"
    )
    parser.add_argument(
        "--center",
        help="the center of the scoop, if ligand is not available, either a "
        "string or a file with the coordinates",
        default="0.0 0.0 0.0",
    )
    parser.add_argument(
        "--innercut",
        type=float,
        help="maximum distance from ligand defining inner region of the scoop",
        default=16.0,
    )
    parser.add_argument(
        "--outercut",
        type=float,
        help="maximum distance from ligand defining outer region of the scoop",
        default=20.0,
    )
    parser.add_argument(
        "--flexin",
        choices=["sidechain", "flexible", "rigid"],
        help="the flexibility of the inner region",
        default="flexible",
    )
    parser.add_argument(
        "--flexout",
        choices=["sidechain", "flexible", "rigid"],
        help="the flexibility of the inner region",
        default="sidechain",
    )
    parser.add_argument(
        "--terminal",
        choices=["keep", "doublekeep", "neutralize"],
        help="controls of to deal with charged terminal",
        default="neutralize",
    )
    parser.add_argument(
        "--excluded",
        nargs="+",
        type=int,
        help="a list of indices for residues to be excluded from scoops",
        default=[],
    )
    parser.add_argument(
        "--added",
        nargs="+",
        type=int,
        help="a list of indices for residues to be included in outer scoops",
        default=[],
    )
    parser.add_argument(
        "--scooplimit",
        help="the minimum difference between number of residues in protein"
        " and scoop for scoop to be retained",
        default=10,
    )
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("scoop_py.log")

    if args.ligand is None:
        ligand = args.center
    else:
        ligand = simulationobjects.PDBFile(filename=args.ligand)
    protein = simulationobjects.PDBFile(filename=args.protein)
    protein = scoop(
        protein,
        ligand,
        args.innercut,
        args.outercut,
        args.flexin,
        args.flexout,
        args.terminal,
        args.excluded,
        args.added,
        args.scooplimit,
    )
    protein.write(args.out, renumber=True)
