#!/usr/bin/env python2.7
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

import simulationobjects as sim
import optimise_gcmcbox

logger = logging.getLogger('protoms')


def make_gcmcbox(pdb, filename, padding=2.0, heavy=False):
    """
    Make a GCMC/JAWS-1 simulation box around a PDB-structure

    Parameters
    ----------
    pdb : PDBFile object
      the PDB structure
    filename : string
      the name of the output file
    padding : float, optional
      the amount of extra space around the ligand to add
    heavy : boolean, optional
      decides to ignore hydrogens, when set to true

    Returns
    -------
    volume : float
      volume of the gcmc region
    """

    logger.debug("Running make_gcmcbox with arguments: ")
    logger.debug("\tpdb      = %s" % pdb)
    logger.debug("\tfilename = %s" % filename)
    logger.debug("\tpadding  = %d" % padding)
    logger.debug("This will make a simulation box for GCMC/JAWS-1")
    print 'Making GCMC box'
    # Create a box around the solute and pad it with two Angstromgs
    box = pdb.getBox(heavy=heavy)
    box["origin"] = box["origin"] - padding
    box["len"] = box["len"] + 2.0 * padding

    volume = volume_box(box["len"])

    # Save it to disc
    sim.write_box(filename, box)
    return volume


def make_gcmcsphere(pdb, filename, padding=2.0, heavy=False):
    """
    Make a GCMC/JAWS-1 simulation sphere around a PDB-structure

    Parameters
    ----------
    pdb : PDBFile object
      the PDB structure
    filename : string
      the name of the output file
    padding : float, optional
      the amount of extra space around the ligand to add

    Returns
    -------
    volume : float
      volume of the gcmc region
    """

    logger.debug("Running make_gcmsphere with arguments: ")
    logger.debug("\tpdb      = %s" % pdb)
    logger.debug("\tfilename = %s" % filename)
    logger.debug("\tpadding  = %d" % padding)
    logger.debug("This will make a simulation box for GCMC/JAWS-1")
    print 'Making GCMC sphere'

    # finds the center of the GCMC sphere
    center, radius = pdb.getSphere(heavy=heavy)
    radius = radius + padding

    sim.write_sphere(filename, center, radius)

    return volume_sphere(radius)


def make_gcmc_region(pdb, box_args, filename, sphere,
                     padding=2.0, heavy=False):
    if sphere:
        volume = make_gcmcsphere(pdb, filename, padding, heavy)
    else:
        if box_args is not None:
            if len(box_args) == 3:
                box = {"len": np.array([args.padding * 2] * 3)}
            elif len(box_args) == 6:
                box = {"len": np.array(map(float, box_args[3:]))}
            else:
                raise sim.SetupError(
                    "Error with 'box' arguement. Please specify either three"
                    " arguments for the centre of the box, or six arguements "
                    "for the centre of the box AND the lengths of the sides.")

            box["center"] = np.array(map(float, box_args[:3]))
            sim.write_box(filename, box)
            volume = volume_box(box["len"])
        else:
            volume = make_gcmcbox(pdb, args.out, args.padding, args.heavy)

    print "Volume of GCMC region:", np.round(volume, 2)
    print "Bequil:", np.round(get_bequil(volume), 2)
    return volume


def volume_box(boxlen):
    return boxlen[0] * boxlen[1] * boxlen[2]


def volume_sphere(radius):
    return (4. / 3.) * np.pi * (radius**3.)


def get_bequil(volume):
    betamu = -10.47
    return betamu + np.log(volume / 30.0)


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to make a PDB-file with box coordinates"
                    " covering a solute molecules")
    parser.add_argument(
        '-s', '--solute',
        help="the name of the PDB-file containing the solute.")
    parser.add_argument(
        '-p', '--padding', type=float,
        help="the padding of box or radius of sphere in A,default=2",
        default=2.0)
    parser.add_argument(
        '-o', '--out',
        help="the name of the box PDB-file",
        default=None)
    parser.add_argument(
        '--sphere',
        help="If flag given, a gcmc sphere will be defined",
        action='store_true')
    parser.add_argument(
        '--heavy',
        help="If given, the GCMC region is built around heavy atoms only",
        action='store_true')
    parser.add_argument(
        '-b', '--box', nargs='+',
        help="Either the centre of the box (x,y,z), or the centre of box AND "
             "length (x,y,z,x,y,z). If the centre is specified and the length "
             "isn't, twice the 'padding' will be the lengths of a cubic box.",
        default=None)
    parser.add_argument(
        '--optimise',
        help="When true, systematically rotates the frame of reference of "
             "the ligand to minimise the box volume. Not relevant when "
             "using a GCMC sphere. Makse sure to supply any other input "
             "ligands or protein(s) to the -f flag so that the same "
             "rotation can be applied.",
        action='store_true')
    parser.add_argument(
      '-f', '--otherfiles', nargs='+',
      help="Other input files (ligands, proteins, etc.) which must be rotated "
           "in the same way to minimise the box size", default=None)
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = sim.setup_logger("make_gcmcbox_py.log")

    if args.solute is None and args.box is None:
        print "Nothing to do! Exiting."

    if args.out is None:
        if args.sphere:
            filename = "gcmc_sphere.pdb"
        else:
            filename = "gcmc_box.pdb"

    try:
        solute_pdbobj = sim.PDBFile()
        solute_pdbobj.read(args.solute)
    except IOError:
        solute_pdbobj = None

    if args.optimise:
        if args.otherfiles is None:
            print('\nNo protein or other input data has been supplied!'
                  'Check that this is not a mistake...\n')
        optimise_gcmcbox.optimise_system(
          ligand_file=args.solute,
          other_files=args.otherfiles,
          heavy=args.heavy)

    volume = make_gcmc_region(solute_pdbobj, args.box, filename,
                              args.sphere, args.padding, args.heavy)


    # if args.sphere:
    #     # volume = make_gcmcsphere(
    #     #     pdbobj, 'gcmc_sphere.pdb', args.padding, args.heavy)
    #     pass
    # else:
    #     if args.box is None:
    #         # pdbobj = sim.PDBFile()
    #         # pdbobj.read(args.solute)
    #         # volume = make_gcmcbox(pdbobj, args.out, args.padding, args.heavy)
    #         pass
    #     elif len(args.box) == 3:
    #         pass
    #         # box = {"center": np.array([float(args.box[0]),
    #         #                            float(args.box[1]),
    #         #                            float(args.box[2])]),
    #         #        "len": np.array([args.padding * 2] * 3)}
    #         # volume = volume_box(box["len"])
    #         # # print_bequil(volume)
    #         # sim.write_box(args.out, box)
    #     elif len(args.box) == 6:
    #         # box = {"center": np.array([float(args.box[0]),
    #         #                            float(args.box[1]),
    #         #                            float(args.box[2])]),
    #         #        "len": np.array([float(args.box[3]),
    #         #                         float(args.box[4]),
    #         #                         float(args.box[5])])}
    #         # volume = volume_box(box["len"])
    #         # print_bequil(volume)
    #         # sim.write_box(args.out, box)
    #         pass
    #     else:
    #         pass
    #         # print ("\nError with 'box' arguement. Please specify either three "
    #         #        "arguements for the centre of the box, or six arguements "
    #         #        "for the centre of the box AND the lengths of the sides.\n")
