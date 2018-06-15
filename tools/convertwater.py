# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines to convert water molecules to water models from a PDB file.

This module defines a single public function
convertwater

Can be executed from the command line as a stand-alone program
"""
from __future__ import print_function
import logging

import numpy as np

from protomslib import simulationobjects
from protomslib.prepare import convertwater

logger = logging.getLogger('protoms')


# --------------------------------------------------
#  Functions to be used for structural alignment of
#  water molecules. Called by rottranstemplate.
# --------------------------------------------------


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to convert water molecules - with or without"
                    " hydrogens - in a pdb file to simulation models, such"
                    " as tip4p. Currently ignores original hydrogen positions."
    )
    parser.add_argument(
        '-p',
        '--pdb',
        help="the PDF-file containing the waters to be transformed")
    parser.add_argument(
        '-o',
        '--out',
        help="the output PDB-file",
        default="convertedwater.pdb")
    parser.add_argument(
        '-m', '--model', help="the water model,default=tip4p", default="tip4p")
    parser.add_argument(
        '-i',
        '--ignoreh',
        action='store_true',
        help="whether to ignore hydrogens in input water. If no hydrogens"
             " are present, waters are randomly orientated. default=No",
        default=False)
    parser.add_argument(
        '-n',
        '--resname',
        help="the residue name that will be applied to the water molecules. "
             "When it is not specified, it is chosen based on the water model",
        default=None)
    parser.add_argument(
        '--setupseed',
        help="optional random number seed for generation of water coordinates",
        default=None,
        type=int)
    return parser


#
# -------------------------------------
# If this is run from the command-line
# -------------------------------------
#

if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("convertwater_py.log")
    if args.setupseed is not None:
        logger.debug("Setup seed = %d" % args.setupseed)
    np.random.seed(args.setupseed)

    pdb_in = simulationobjects.PDBFile(filename=args.pdb)
    pdb_out = convertwater(
        pdb_in, args.model, ignorH=args.ignoreh, watresname=args.resname
    )
    pdb_out.write(args.out)
