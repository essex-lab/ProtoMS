# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Classes and routines to make ProtoMS command files

This module defines a single public function
generate_input

and the following public classes:
ProtoMSSimulation
ProteinLigandSimulation
Equilibration
Sampling
DualTopology

Can be executed from the command line as a stand-alone program
"""

import logging
from protomslib import simulationobjects
from protomslib.command import generate_input

logger = logging.getLogger('protoms')


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to create a ProtoMS command file")
    parser.add_argument(
        '-s',
        '--simulation',
        choices=[
            "sampling", "equilibration", "dualtopology", "singletopology",
            "gcap_single","gcap_dual","gcmc", "jaws1", "jaws2"
        ],
        help="the kind of simulation to setup",
        default="equilibration")
    parser.add_argument(
        '--dovacuum',
        action='store_true',
        help="turn on vacuum simulation for simulation types equilibration"
             " and sampling",
        default=False)
    parser.add_argument('-p', '--protein', help="the name of the protein file")
    parser.add_argument(
        '-l', '--ligands', nargs="+", help="the name of the ligand pdb files")
    parser.add_argument(
        '-t',
        '--templates',
        nargs="+",
        help="the name of ProtoMS template files")
    parser.add_argument(
        '-pw', '--protwater', help="the name of the solvent for protein")
    parser.add_argument(
        '-lw', '--ligwater', help="the name of the solvent for ligand")
    parser.add_argument(
        '-o',
        '--out',
        help="the prefix of the name of the command file",
        default="run")
    parser.add_argument(
        '--outfolder', help="the ProtoMS output folder", default="out")
    parser.add_argument(
        '--gaff',
        help="the version of GAFF to use for ligand",
        default="gaff16")
    parser.add_argument(
        '--lambdas',
        nargs="+",
        type=float,
        help="the lambda values or the number of lambdas",
        default=[16])
    parser.add_argument(
        '--adams',
        nargs="+",
        type=float,
        help="the Adam/B values for the GCMC",
        default=0)
    parser.add_argument(
        '--adamsrange',
        nargs="+",
        type=float,
        help="the upper and lower Adam/B values for the GCMC and, optionally,"
             " the number of values desired (default value every 1.0), e.g. -1"
             " -16 gives all integers between and including -1 and -16",
        default=None)
    parser.add_argument(
        '--jawsbias',
        nargs="+",
        type=float,
        help="the bias for the JAWS-2",
        default=0)
    parser.add_argument(
        '--gcmcwater', help="a pdb file with a box of water to do GCMC on")
    parser.add_argument(
        '--gcmcbox', help="a pdb file with box dimensions for the GCMC box")
    parser.add_argument(
        '--watmodel',
        help="the name of the water model. Default = tip4p",
        choices=['tip3p', 'tip4p'],
        default='tip4p')
    parser.add_argument(
        '--nequil',
        type=float,
        help="the number of equilibration steps",
        default=5E6)
    parser.add_argument(
        '--nprod',
        type=float,
        help="the number of production steps",
        default=40E6)
    parser.add_argument(
        '--dumpfreq',
        type=float,
        help="the output dump frequency",
        default=1E5)
    parser.add_argument(
        '--absolute',
        action='store_true',
        help="whether an absolute free energy calculation is to be run."
             " Default=False",
        default=False)
    parser.add_argument(
        '--ranseed',
        help="the value of the random seed you wish to simulate with. "
             "If None, then a seed is randomly generated. Default=None",
        default=None)
    parser.add_argument(
        '--tune', action='store_true', help=argparse.SUPPRESS, default=False)
    parser.add_argument(
        '--softcore', type=str, default='all',
        choices=('auto', 'all', 'none', 'manual'),
        help="determine which atoms to apply softcore potentials to. If 'all' "
             "softcores are applied to all atoms of both solutes. If 'none' "
             "softcores are not applied to any atoms. If 'auto', softcores are"
             " applied to atoms based on matching coordinates between ligand "
             "structures. The selected softcore atoms can be amended using the"
             " --spec-softcore flag. If 'manual' only those atoms specified by"
             " the --spec-softcore flag are softcore.")
    parser.add_argument(
        '--spec-softcore', type=str,
        help='Specify atoms to add or remove from softcore selections. Can be '
             'up to two, space separated, strings of the form "N:AT1,AT2,-AT3"'
             '. N should be either "1" or "2" indicating the corresponding '
             'ligand. The comma separated list of atom names are added to the'
             ' softcore selection. A preceding dash for an atom name specifies'
             ' it should be removed from the softcore selection. The special '
             'value "auto" indictates that automatic softcore assignments '
             'should be accepted without amendment.')
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("generate_input_py.log")

    free_cmd, bnd_cmd, gas_cmd = generate_input(
        args.protein, args.ligands, args.templates, args.protwater,
        args.ligwater, args.ranseed, args)

    # protoMS cannot handle cmd files containing upper case letters
    args.out = args.out.lower()
    if free_cmd is not None:
        free_cmd.writeCommandFile(args.out + "_free.cmd")
    if bnd_cmd is not None:
        if args.simulation == "gcmc":
            bnd_cmd.writeCommandFile(args.out + "_gcmc.cmd")
        elif args.simulation in ["jaws1", "jaws2"]:
            bnd_cmd.writeCommandFile(args.out + "_jaws.cmd")
        else:
            bnd_cmd.writeCommandFile(args.out + "_bnd.cmd")
    if gas_cmd is not None:
        gas_cmd.writeCommandFile(args.out + "_gas.cmd")

