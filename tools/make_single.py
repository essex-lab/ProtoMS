# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Routines to make ProtoMS single topology template

This module defines these public function
make_single
write_map
summarize

Can be executed from the command line as a stand-alone program
"""

import logging

from protomslib import simulationobjects as sim
from protomslib.templates import make_single, write_map, summarize_single

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to setup template files for single-toplogy "
        "perturbations semi-automatically"
    )
    parser.add_argument("-t0", "--tem0", help="Template file for V0")
    parser.add_argument("-t1", "--tem1", help="Template file for V1")
    parser.add_argument("-p0", "--pdb0", help="PDB-file for V0")
    parser.add_argument("-p1", "--pdb1", help="PDB-file for V1")
    parser.add_argument(
        "-m", "--map", help="the correspondance map from V0 to V1"
    )
    parser.add_argument(
        "-o", "--out", default="single", help="prefix of the output file"
    )
    parser.add_argument(
        "--gaff",
        default="gaff16",
        help="the version of GAFF to use for ligand",
    )
    return parser


#
# If this is run from the command-line
#
if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = sim.setup_logger("make_single_py.log")

    if None in (args.tem0, args.tem1, args.pdb0, args.pdb1):
        sim.SetupError(
            "Not all four necessary input files given. Cannot setup single-topology"
        )

    eletem, vdwtem, combtem, cmap = make_single(
        args.tem0,
        args.tem1,
        args.pdb0,
        args.pdb1,
        args.map,
        gaffversion=args.gaff,
    )
    eletem.write(args.out + "_ele.tem")
    vdwtem.write(args.out + "_vdw.tem")
    combtem.write(args.out + "_comb.tem")
    if args.map is None:
        write_map(cmap, args.out + "_cmap.dat")

    # Write out a summary
    summarize_single(eletem, vdwtem, logger.info)
