# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
""" Routines to convert atom names in a PDB-file

This module defines two public functions:
read_convfile
pdb2pms

Can be executed from the command line as a stand-alone program
"""

import logging
from protomslib import simulationobjects
from protomslib.prepare import pdb2pms

logger = logging.getLogger("protoms")


def get_arg_parser():

    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program convert atom names in a protein pdb-file to"
        " ProtoMS style"
    )
    parser.add_argument("-p", "--protein", help="the protein PDB-file")
    parser.add_argument(
        "-o", "--out", help="the output PDB-file", default="protein_pms.pdb"
    )
    parser.add_argument(
        "-s",
        "--style",
        help="the style of the input PDB-file",
        default="amber",
    )
    parser.add_argument(
        "-c",
        "--conversionfile",
        help="the name of the file with conversion rules",
        default="atomnamesmap.dat",
    )
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("convertatomnames_py.log")

    protein = simulationobjects.PDBFile(filename=args.protein)
    protein_out = pdb2pms(protein, args.style, args.conversionfile)
    protein_out.write(args.out)
