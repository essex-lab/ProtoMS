# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines to solvate a structure

This module defines a single public function:
solvate

Can be executed from the command line as a stand-alone program
"""
import numpy as np
import six
from protomslib import simulationobjects
from protomslib.solvate import solvate


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    disclaimer = """
    if -b or -s are not supplied on the command-line, the program will ask for
    them.\n
    -c can be either 'cent' or a string containing 1, 2 or 3 numbers. If 1
    number is given it will be used as center of the droplet in x, y, and z.
    If 2 numbers are given this is interpreted as an atom range, such that the
    droplet will be centered on the indicated atoms, and if 3 numbers are
    given this is directly taken as the center of droplet\n
    Example usages:
      solvate.py -b ${PROTOMSHOME}/tools/sbox1.pdb -s solute.pdb
         (will solvate 'solute.pdb' in a box that extends at least 10 A from
          the solute)
      solvate.py -b ${PROTOMSHOME}/tools/sbox1.pdb -s protein.pdb -g droplet -r 25.0
         (will solvate 'protein.pdb' in a 25 A droplet centered on
          all coordinates)
  """
    parser = argparse.ArgumentParser(
        description="Program to solvate a solute molecule in either a box or "
                    "a droplet",
        epilog=disclaimer,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-b',
        '--box',
        help="a PDB-file containing a pre-equilibrated box of water molcules",
        default="")
    parser.add_argument(
        '-s',
        '--solute',
        help="a PDB-file containing the solute molecule",
        default=None)
    parser.add_argument(
        '-pr',
        '--protein',
        help="a PDB-file containing the protein molecule",
        default=None)
    parser.add_argument(
        '-o',
        '--out',
        help="the name of the output PDB-file containing the added water, "
             "default solvent_box.pdb",
        default="solvent_box.pdb")
    parser.add_argument(
        '-g',
        '--geometry',
        choices=["box", "droplet", "flood"],
        help="the geometry of the added water, should be either 'box',"
             " 'droplet' or 'flood'",
        default="box")
    parser.add_argument(
        '-p',
        '--padding',
        type=float,
        help="the minimum distance between the solute and the box edge,"
             " default=10 A",
        default=10.0)
    parser.add_argument(
        '-r',
        '--radius',
        type=float,
        help="the radius of the droplet, default=30A",
        default=30.0)
    parser.add_argument(
        '-c',
        '--center',
        help="definition of center, default='cent'",
        default="cent")
    parser.add_argument(
        '-n',
        '--names',
        choices=["Amber", "ProtoMS"],
        help="the naming convention, should be either Amber or ProtoMS",
        default="ProtoMS")
    parser.add_argument(
        '--offset',
        type=float,
        help="the offset to be added to vdW radii of the atoms to avoid"
             " overfilling cavities with water.",
        default=0.89)
    parser.add_argument(
        '--setupseed',
        type=int,
        help="optional random number seed for generation of water"
             " coordinates..",
        default=None)
    return parser


#
# -------------------------------------
# If this is run from the command-line
# -------------------------------------
#
if __name__ == '__main__':

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("solvate_py.log")
    if args.setupseed is not None:
        logger.debug("Setup seed = %d" % args.setupseed)
    np.random.seed(args.setupseed)

    # Ask for input that is absolutely necessary and
    # that do not have any defaults
    if args.solute is None and args.protein is None:
        print("You haven't entered a protein or a solute to solvate!\n")
        print("Please enter a solute now (Return to skip):")
        args.solute = six.moves.input()
        print("Please enter a protein now (Return to skip):")
        args.protein = six.moves.input()
        if args.solute == "" and args.protein == "":
            print(
                "You still haven't entered a protein or a solute to solvate!",
                end='\n\n')
            quit()
    if args.box == "":
        print("Enter the filename of a pre-equilibrated water box: ", end="")
        args.box = six.moves.input()
    if args.solute == "":
        args.solute = None
    if args.protein == "":
        args.protein = None

    boxpdb = solvate(args.box, args.solute, args.protein, args.geometry,
                     args.padding, args.radius, args.center, args.names,
                     args.offset)
    boxpdb.write(args.out)
