#!/usr/bin/env python
# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Program to setup, run and analyse a ProtoMS simulation
"""

from __future__ import print_function
import os
import sys
import time
import numpy as np
import six
from protomslib import simulationobjects
from protomslib.prepare import _get_prefix, _load_ligand_pdb, _prep_ligand
from protomslib.prepare import _prep_gcmc, _prep_protein, _prep_singletopology
from protomslib.prepare import _merge_templates, _prep_jaws2, make_dummy
from protomslib.command import generate_input
from protomslib.utils import _cleanup


# Logger is used globally
logger = simulationobjects.setup_logger("protoms_py.log")


def _wizard(settings):

    print("\nBecause you initiated protoms.py without any arguments,"
          "a wizard will walk you through the setup\n")

    print("What kind of simulation would you like to setup?")
    print("\t1) Equilibration")
    print("\t2) Sampling")
    print("\t3) Dual topology free energy")
    print("\t4) Single topology free energy")
    print("\t5) Grand-canonical Monte Carlo (GCMC)")
    print("\t6) JAWS stage 1 (Just Add Waters)")
    print("\t7) JAWS stage 2 (Just Add Waters)")
    print(">", end="")
    valid = ["", "1", "2", "3", "4", "5", "6", "7"]
    instr = six.moves.input()
    if instr in ["jon", "dev"]:
        import matplotlib.pylab as plt
        img = np.load(simulationobjects.standard_filename(".ee.npz",
                                                          "tools"))[instr]
        plt.imshow(img)
        plt.show()
        return
    while instr not in valid:
        print("Please type a number between 1 and 7!")
        print(">", end="")
        instr = six.moves.input()
    if instr == "":
        return
    val = int(instr)
    vals = [
        "equilibration", "sampling", "dualtopology", "singletopology", "gcmc",
        "jaws1", "jaws2"
    ]
    settings.simulation = vals[val - 1]

    print("\nDo you have a protein that you would like to setup and simulate?")
    print("\tPlease enter the protein name or a PDB filename.\n\tPress enter"
          "if you don't have a protein.\n>", end="")
    instr = six.moves.input()
    if instr != "":
        settings.protein = instr

    print("\nDo you have a ligand that you would like to setup and simulate?")
    if settings.simulation not in ["dualtopology", "singletopology"]:
        print("\tPlease enter the ligand name or a PDB filename.")
        print("\tPress enter if you don't have a ligand.\n>", end="")
        instr = six.moves.input()
        if instr != "":
            args.ligand = [instr]
    else:
        while instr == "":
            print("For the type of simulation you have selected, two ligands "
                  " need to be given. \n\tPlease enter a ligand name or a PDB"
                  " filename\n>", end="")
            instr = six.moves.input()
        args.ligand = [instr]
        print("\tPlease enter the second ligand name or a PDB filename.")
        instr = ""
        while instr == "":
            print(">", end="")
            instr = six.moves.input()
            if instr == "":
                print(
                    "\tPlease note that absolute free energies cannot "
                    "currently be set up via the wizard.\n\tPlease see "
                    "protoms.py --fullhelp for details of how to set these up!"
                )
                print(
                    "\tPlease enter the second ligand name or a PDB filename.")
        args.ligand.append(instr)

    print("\nDo you have a co-factor or another solute that you would like "
          "to setup and simulate?")
    print("\tPlease enter the solute name or a PDB filename.")
    print("\tPress enter if you don't have a solute or when you have "
          "specified all solutes.")
    print(">", end="")
    instr = six.moves.input()
    if instr != "":
        if args.ligand is None or len(args.ligand) == 0:
            args.ligand = [instr]
        else:
            args.ligand.append(instr)
    while instr != "":
        print(">", end="")
        instr = six.moves.input()
        if instr != "":
            args.ligand.append(instr)


if __name__ == "__main__":

    # Setup a parser of the command-line arguments
    parser = simulationobjects.MyArgumentParser(
        description="Program setup and run a ProtoMS simulations")
    parser.add_argument(
        '-s', '--simulation', default="none",
        choices=["none", "equilibration", "sampling", "dualtopology",
                 "singletopology", "gcmc", "jaws1", "jaws2"],
        help="the kind of simulation to setup")
    parser.add_argument(
        '-f', '--folders', nargs="+", default=["."],
        help="folders to search for files ")
    parser.add_argument('-p', '--protein', help="the prefix of the protein")
    parser.add_argument(
        '-l', '--ligand', nargs="+", help="the prefix of the ligand(s)")
    parser.add_argument(
        '-w', '--water', default="water",
        help="the prefix of the water/solvent")
    parser.add_argument(
        '-c', '--cmdfile', default="run",
        help="the prefix of the command file")
    parser.add_argument(
        '-sc', '--scoop', help="the name of your protein scoop")
    parser.add_argument(
        '-t', '--template', nargs="+",
        help="the template files for your ligands")
    parser.add_argument(
        '-r', '--repeats', default="",
        help="the number of repeats to be run (if more than 1) or a name"
             " for your repeat")
    # General control variables
    cntrlgroup = parser.add_argument_group("General control variables")
    cntrlgroup.add_argument(
        '--outfolder', help="the ProtoMS output folder", default="")
    cntrlgroup.add_argument(
        '--atomnames', help="a file with atom name conversions")
    cntrlgroup.add_argument(
        '--watmodel', choices=['tip3p', 'tip4p'], default='tip4p',
        help="the name of the water model. Default = tip4p")
    cntrlgroup.add_argument(
        '--waterbox', help="a file with pre-equilibrated water molecules")
    cntrlgroup.add_argument(
        '--setupseed', type=int, default=None,
        help="optional seed for random number generators in setup")
    # Ligand setup variables
    liggroup = parser.add_argument_group("Ligand setup variables")
    liggroup.add_argument(
        '--charge', nargs="+", type=float,
        help="the net charge of each ligand")
    liggroup.add_argument(
        '--singlemap', help="the correspondance map for single-topology")
    liggroup.add_argument(
        '--gaff', default="gaff16",
        help="the version of GAFF to use for ligand")
    # Protein setup variables
    protgroup = parser.add_argument_group("Protein setup variables")
    protgroup.add_argument(
        '--center', default=None,
        help="the center of the scoop, if ligand is not available, either "
             "a string or a file with the coordinates")
    protgroup.add_argument(
        '--innercut', type=float, default=16.0,
        help="maximum distance from ligand defining inner region of the scoop")
    protgroup.add_argument(
        '--outercut', type=float, default=20.0,
        help="maximum distance from ligand defining outer region of the scoop")
    protgroup.add_argument(
        '--flexin', default="flexible",
        choices=['sidechain', 'flexible', 'rigid'],
        help="the flexibility of the inner region")
    protgroup.add_argument(
        '--flexout',
        choices=['sidechain', 'flexible', 'rigid'],
        help="the flexibility of the outer region",
        default="sidechain")
    protgroup.add_argument(
        '--scooplimit',
        help="the minimum difference between number of residues in protein "
             "and scoop for scoop to be retained",
        default=10)
    protgroup.add_argument(
        '--capradius', type=float, default=30.0,
        help="the radius of the droplet around the protein")
    # Simulation parameters
    simgroup = parser.add_argument_group("Simulatiom parameters")
    simgroup.add_argument(
        '--lambdas',
        nargs="+",
        type=float,
        help="the lambda values or the number of lambdas",
        default=[16])
    simgroup.add_argument(
        '--adams',
        nargs="+",
        type=float,
        help="the Adam/B values for the GCMC",
        default=0)
    simgroup.add_argument(
        '--adamsrange', nargs="+", type=float, default=None,
        help="the upper and lower Adam/B values for the GCMC and, optionally"
             ", the number of values desired (default value every 1.0), e.g."
             " -1 -16 gives all integers between and including -1 and -16")
    simgroup.add_argument(
        '--gcmcwater',
        help="a pdb file with a box of water to do GCMC on or an integer "
             "corresponding to the number of water molecules to add"
    )
    simgroup.add_argument(
        '--gcmcbox', nargs="+",
        help="a pdb file with box dimensions for the GCMC box, or a list "
             "of origin(x,y,z) and length(x,y,z) coordinates"
    )
    simgroup.add_argument(
        '--jawsbias', nargs="+", type=float, default=[6.5],
        help="the bias in JAWS-2")
    simgroup.add_argument(
        '--nequil', type=float, default=5E6,
        help="the number of equilibration steps")
    simgroup.add_argument(
        '--nprod', type=float, default=40E6,
        help="the number of production steps")
    simgroup.add_argument(
        '--dumpfreq', type=float, default=1E5,
        help="the output dump frequency")
    simgroup.add_argument(
        '--ranseed', default=None,
        help="the value of the random seed you wish to simulate with. "
             "If None, then a seed is randomly generated. Default=None")
    simgroup.add_argument(
        '--absolute', action='store_true', default=False,
        help="whether an absolute free energy calculation is to be run. "
             "Default=False")
    simgroup.add_argument(
        '--dovacuum', action='store_true', default=False,
        help="turn on vacuum simulation for simulation types equilibration"
             " and sampling")
    simgroup.add_argument(
        '--testrun', action='store_true', default=False,
        help="setup a short test run. Default=False")
    simgroup.add_argument(
        '--cleanup', action='store_true', default=False,
        help="Clean up extra files. Default=False")
    simgroup.add_argument(
        '--tune', action='store_true', default=False,
        help='Carry out dihedral tuning simulation')
    simgroup.add_argument(
        '--softcore', type=str, default='all',
        choices=('auto', 'all', 'none', 'manual'),
        help="determine which atoms to apply softcore potentials to.\n "
             "'all'=softcores applied to all atoms of both solutes, "
             "'none'=softcores not applied to any atoms\n "
             "'mixed'=softcores will be applied only to non matching "
             "atoms within ligand structures")
    parser.add_argument(
      '--spec-softcore', type=str,
      help='Specify atoms to add or remove from softcore selections. Can be '
           'up to two, space separated, strings of the form "N:AT1,AT2,-AT3". '
           'N should be either "1" or "2" indicating the corresponding ligand.'
           ' The comma separated list of atom names are added to the softcore '
           'selection. A preceding dash for an atom name specifies it should '
           'be removed from the softcore selection.')
    args = parser.parse_args()

    print(r"""
            __|\
         .-'    '-.
        / .--, _ a L
      .J (  '-' "'--'
     '-'-.)  .~~~~~~~~~~~~~~~~~~~~.
             |                    |     __
             |    Welcome to      | ,.-'e ''-'7
             |    protoms.py      |  '--.    (
             |                    |      ),   \
             '~~~~~~~~~~~~~~~~~~~~'      ` )  :
                                      ,__.'_.'
                                      '-, (
                                        '--'  """)
    print()

    # Setup the logger - logger is global
    logger.debug(
        "Running protoms.py at %s" % time.strftime("%d/%m/%Y - %H:%M:%S"))
    logger.debug("Command line arguments = %s" % " ".join(sys.argv[1:]))
    logger.debug("Settings = %s" % args)

    np.random.seed(args.setupseed)

    # Adds current folder to the folders
    args.folders.append(".")

    # If we run without any command-line arguments, initiate the wizard
    if len(sys.argv) == 1:
        _wizard(args)

    if args.protein is None and args.ligand is None and args.scoop is None:
        print("Nothing to do, so exit!")
        quit()

    # Set $PROTOMSHOME
    if os.getenv("PROTOMSHOME") is None:
        string = os.path.dirname(os.path.abspath(__file__))
        logger.info("Setting PROTOMSHOME to %s" % string)
        os.environ[
            "PROTOMSHOME"] = string  # This does not change the original shell

    # Try to find a default water box
    if args.waterbox is None:
        args.waterbox = simulationobjects.standard_filename(
            "wbox_" + args.watmodel.lower() + ".pdb", "data")
    if not os.path.isfile(args.waterbox):
        msg = "Could not find file (%s) with pre-equilibrated waters" % \
            args.waterbox
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    # Setup list of files to be stored away
    tarlist = []

    # Prepare each given ligand
    ligand_files = {
    }  # This will be filled with a dictionary of filenames for each ligand
    ligands = []  # This will be a list of ligands
    ligpdbs = None  # This will be a list of ligand pdb-files
    ligtems = None  # This will be a list of ligand template files
    ligobjs = None  # This will hold a merged pdb object of all ligand pdbs
    ligand_water = None  # This will hold the filename of the free-leg waterbox
    if args.ligand is not None:
        # Read in each ligand pdb file and create a pdb object
        for l in args.ligand:
            prefix = _get_prefix(l)
            ligands.append(prefix)
            ligand_files[prefix] = {}
            ligand_files[prefix]["pdb"], ligand_files[prefix][
                "obj"] = _load_ligand_pdb(prefix, args.folders)

        # Make a dummy PDB structure
        if args.simulation == "dualtopology" and len(ligands) < 2:
            if args.absolute:
                ligands.insert(1, "*dummy")
                prefix0 = ligands[0]
                dummy_name = os.path.basename(
                    _get_prefix(prefix0)) + "_dummy.pdb"
                ligand_files["*dummy"] = {}
                ligand_files["*dummy"]["pdb"] = dummy_name
                ligand_files["*dummy"]["obj"] = make_dummy(
                    ligand_files[prefix0]["obj"])
                ligand_files["*dummy"]["obj"].write(
                    ligand_files["*dummy"]["pdb"])
                ligand_files["*dummy"][
                    "tem"] = simulationobjects.standard_filename(
                        "dummy.tem", "data")
                logger.info("")
                logger.info("Creating dummy PDB-file for ligand: %s" %
                            ligand_files["*dummy"]["pdb"])
            else:
                msg = "You are trying to run a dual topology simulation " \
                      "but with only one ligand! (%s.pdb)\nIf you meant " \
                      "to run an absolute free energy calculation, rerun " \
                      "with the --absolute option." % ligands[0]
                logger.error(msg)
                raise simulationobjects.SetupError(msg)

        # Create merged pdb objects
        if len(ligands) >= 2:
            ligobj12 = simulationobjects.merge_pdbs(
                ligand_files[l]["obj"] for l in ligands[:2])
        else:
            ligobj12 = ligand_files[ligands[0]]["obj"]
        if len(ligands) > 1:
            ligobjs = simulationobjects.merge_pdbs(
                ligand_files[l]["obj"] for l in ligands)
        else:
            ligobjs = ligand_files[ligands[0]]["obj"]

        # Now do the preparations
        for i, l in enumerate(args.ligand):
            if l[0] == "*":
                continue  # Skip ligands created in the script, i.e. the dummy
            if args.charge is not None and i < len(args.charge):
                charge = args.charge[i]
            else:
                charge = 0
            if i > 1:
                ligobj12 = None
            prefix = _get_prefix(l)
            _prep_ligand(ligand_files[prefix], i == 0, charge, ligobj12,
                         args.folders, tarlist, args)

        ligpdbs = [ligand_files[l]["pdb"] for l in ligands]
        ligtems = [ligand_files[l]["tem"] for l in ligands]
        ligand_water = ligand_files[ligands[0]]["wat"]

        # Here we need to make single topology templates, if requested
        if args.simulation == "singletopology":
            ligtems, ligtems2, ligtems3 = _prep_singletopology(
                ligpdbs, ligtems, tarlist, args)

        # Here we will merge ligand template files if there is more than one
        if len(ligtems) > 1:
            logger.info("")
            ligtems = _merge_templates(ligtems, tarlist)
            if args.simulation == "singletopology":
                ligtems2 = _merge_templates(ligtems2, tarlist)
                ligtems3 = _merge_templates(ligtems3, tarlist)

    # Prepare the protein
    protein_file = None
    water_file = None
    if args.protein is not None or args.scoop is not None:
        protein_file, water_file = _prep_protein(
            args.protein, ligobjs, args.water, args.folders, tarlist, args)

    # Extra preparation for GCMC or JAWS-1
    if args.simulation in ["gcmc", "jaws1"]:
        if water_file is None:
            msg = "GCMC and JAWS1 not supported without protein or scoop"
            logger.error(msg)
            raise simulationobjects.SetupError(msg)
        args.gcmcwater, water_file = _prep_gcmc(ligands, ligand_files,
                                                water_file, tarlist, args)

    # Extra preparation for JAWS-2
    if args.simulation == "jaws2":
        single_wat, other_wat, water_file = _prep_jaws2(
            water_file, tarlist, args)

    # Check of test run
    if args.testrun:
        if args.nequil == 5E6:
            args.nequil = 0
        if args.nprod == 40E6:
            args.nprod = 4000
        if args.dumpfreq == 1E5:
            args.dumpfreq = 10
        if len(args.lambdas) == 1 and args.lambdas[0] == 16:
            args.lambdas = [4]

    # Create ProtoMS command files
    ranseed = args.ranseed
    if args.simulation == "singletopology":
        postfix = ["_ele", "_vdw", "_comb"]
    elif args.simulation == "jaws2":
        postfix = ["_jaws2-w%d" % (i + 1) for i in range(len(single_wat.pdbs))]
        if args.outfolder == "":
            args.outfolder = "out"
    else:
        postfix = [""]
    if args.repeats.isdigit():
        args.repeats = range(1, int(args.repeats) + 1)
    else:
        args.repeats = [args.repeats.lower()]

    repeats = []
    for post in postfix:
        for repeat in args.repeats:
            repeats.append(str(repeat) + post)

    outfolder = args.outfolder
    if args.outfolder == "":
        if args.simulation == "gcmc":
            outfolder = "out_gcmc"
        elif args.simulation in ["jaws1", "jaws2"]:
            outfolder = "out_jaws"
        else:
            outfolder = "out"

    for repeat in repeats:
        args.outfolder = outfolder + repeat
        if args.simulation not in ["singletopology", "jaws2"] or \
           "_ele" in repeat:

            free_cmd, bnd_cmd, gas_cmd = generate_input(
                protein_file, ligpdbs, ligtems, water_file, ligand_water,
                ranseed, args)
        elif args.simulation == "singletopology" and "_vdw" in repeat:
            free_cmd, bnd_cmd, gas_cmd = generate_input(
                protein_file, ligpdbs, ligtems2, water_file, ligand_water,
                ranseed, args)
        elif args.simulation == "singletopology" and "_comb" in repeat:
            free_cmd, bnd_cmd, gas_cmd = generate_input(
                protein_file, ligpdbs, ligtems3, water_file, ligand_water,
                ranseed, args)
        elif args.simulation == "jaws2":
            idx = int(repeat.split("-")[-1][1:])
            args.gcmcwater = "jaws2_wat%d.pdb" % idx
            jaws2wat = "jaws2_not%d.pdb" % idx
            free_cmd, bnd_cmd, gas_cmd = generate_input(
                protein_file, ligpdbs, ligtems, water_file + " " + jaws2wat,
                ligand_water, ranseed, args)

        # ProtoMS cannot handle command files with lower case letters
        args.cmdfile = args.cmdfile.lower()
        if free_cmd is not None:
            free_cmd.writeCommandFile(args.cmdfile + repeat + "_free.cmd")
        if bnd_cmd is not None:
            if args.simulation == "gcmc":
                bnd_cmd.writeCommandFile(args.cmdfile + repeat + "_gcmc.cmd")
            elif args.simulation in ["jaws1", "jaws2"]:
                bnd_cmd.writeCommandFile(args.cmdfile + repeat + "_jaws.cmd")
            else:
                bnd_cmd.writeCommandFile(args.cmdfile + repeat + "_bnd.cmd")
        if gas_cmd is not None:
            if args.absolute:
                # in this case, gas_cmd contains a cmd file to account for
                # introduction of the harmonic restraint
                gas_cmd.writeCommandFile(args.cmdfile + repeat +
                                         "_bnd_rstr.cmd")
            else:
                gas_cmd.writeCommandFile(args.cmdfile + repeat + "_gas.cmd")

    if args.cleanup:
        _cleanup(tarlist)
