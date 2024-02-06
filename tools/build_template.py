# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines to build a ProtoMS template file

The module defines a two public function
build_template
make_zmat

Can be executed from the command line as a stand-alone program

"""
from __future__ import print_function
import os
import logging

from protomslib import simulationobjects as sim
from protomslib.templates import build_template

logger = logging.getLogger("protoms")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to build a ProtoMS template file"
    )
    parser.add_argument(
        "-p", "--prepi", help="the name of the leap prepi-file"
    )
    parser.add_argument(
        "-o", "--out", help="the name of the template file", default="lig.tem"
    )
    parser.add_argument(
        "-z", "--zmat", help="the name of the zmatrix-file, if it exists"
    )
    parser.add_argument(
        "-f", "--frcmod", help="the name of the frcmod-file, if it exists"
    )
    parser.add_argument(
        "-n", "--name", help="the name of the solute", default="UNK"
    )
    parser.add_argument(
        "-t",
        "--translate",
        help="maxmium size for translation moves in Angstroms",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "-r",
        "--rotate",
        help="maxmium size for rotation moves in degrees",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--alldihs",
        help="sample improper dihedrals",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--gaff",
        help="gaff version to use, gaff14 or gafff16",
        default="gaff16",
    )
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = sim.setup_logger("build_template_py.log")

    tem = build_template(
        temfile=args.out,
        prepifile=args.prepi,
        zmatfile=args.zmat,
        frcmodfile=args.frcmod,
        resname=args.name,
        translate=args.translate,
        rotate=args.rotate,
        alldihs=args.alldihs,
        gaffversion=args.gaff,
    )
    tem.write(args.out)
    if args.zmat is None:
        tem.templates[0].write_zmat(os.path.splitext(args.out)[0] + ".zmat")
