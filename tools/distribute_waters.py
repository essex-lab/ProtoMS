# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Program to generate a set of n molecules contained within given box dimenstions
or to redistribute molecules given in a pdb object or file within the given box dimensions

This program takes the box dimensions as a list of strings convertible to floats,
in the following order :
  origin in x | origin in y | origin in z | length in x | length in y | length in z
And either the number of molecules as a string convertible to integer,
or the pdb, which can be provided as a string with the name of the file, or an object.

It produces a pdb file with the randomly distributed molecules with their oxygen atoms
within the box limits.

Initially thought to generate JAWS or GCMC molecules

Can be executed from the command line as a stand-alone program
"""
import logging
import numpy as np
from protomslib import simulationobjects
from protomslib.gcmc import distribute_particles

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Randomly distribute n molecules within box dimensions"
    )
    parser.add_argument(
        "-b",
        "--box",
        nargs=6,
        help="Dimensions of the box. Six arguments expected: origin (x,y,z) "
        "& length (x,y,z)",
    )
    parser.add_argument(
        "-m",
        "--molecules",
        help="Molecules to distribute in the box. Either the number of waters "
        "or a pdb file containing all of them",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Name of the pdb file to write the molecules to. "
        "Default='ghostmolecules.pdb'",
        default="ghostmolecules.pdb",
    )
    parser.add_argument(
        "--model",
        help="Water model. Used when only the amount of waters is specified. "
        "Options: 't4p','t3p'. Default='t4p'",
        default="t4p",
    )
    parser.add_argument(
        "--resname",
        help="Residue name of the molecules writen to output. Default='WAT'",
        default="WAT",
    )
    parser.add_argument(
        "--number",
        help="Required number of molecules when it differs from the number"
        " of residues in the file.",
        default=None,
    )
    parser.add_argument(
        "--setupseed",
        help="Optional random number seed for generation of water coordinates",
        default=None,
        type=int,
    )
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    logger = simulationobjects.setup_logger("distribute_waters_py.log")
    if args.setupseed is not None:
        logger.debug("Setup seed = %d" % args.setupseed)
    np.random.seed(args.setupseed)

    outobj = distribute_particles(
        args.box, args.molecules, args.model, args.resname, args.number
    )
    outobj.write(filename=args.outfile)
    print("\nMolecules printed in %s" % args.outfile)
