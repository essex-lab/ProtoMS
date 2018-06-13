# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routines that wrap AmberTools programs

This module defines two public functions:
run_antechamber
run_parmchk

Can be executed from the command line as a stand-alone program
"""

import os
import subprocess
import tempfile
import logging
import distutils.spawn
import six

if __name__ != "__main__":
    from . import simulationobjects

logger = logging.getLogger('protoms')


def _run_program(name, command):
    """
    Wrapper for an AmberTools program with default parameters

    run_program ( name, command )

    Parameters
    ----------
    name : string
      the name of the program to run
    command : string
      the command to execute it

    Raises
    ------
    SetupError
      if the program failed to execute properly
    """

    # Create temporary file to write stdout, stderr of the program
    tmpfile, tmpname = tempfile.mkstemp()
    ret_code = subprocess.call(
        command, shell=True, stdout=tmpfile, stderr=tmpfile)
    # Catch some error codes
    if ret_code == 127:
        msg = "Unable to find executable, please make sure this is" \
            " present in your PATH or $AMBERHOME/bin."
        logger.error(msg)
        raise simulationobjects.SetupError(msg)
    if ret_code == 1:
        # Get the error message from the temporary file
        errmsg = "\n".join(line for line in open(tmpname).readlines())
        os.remove(tmpname)
        msg = "%s was not able to run successfully. Please check output. " \
            " Error message was:\n%s" % (name, errmsg)
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    os.remove(tmpname)


def _get_executable_path(name):
    """Robust routine that looks for Amber Tools executables first in
    $AMBERHOME/bin and then in the path generally.

    Parameters
    ----------
    name : string
      the name of the program to find

    Returns
    -------
    string
        the full path of the desired executable

    Raises
    ------
    SetupError
        if the executable could not be found
    """
    try:
        os.environ['AMBERHOME']
    except KeyError:
        raise simulationobjects.SetupError(
            "The environmental variable $AMBERHOME is not set. This is "
            "required to use components of AMBER Tools."
        )

    try:
        exe_path = os.environ['AMBERHOME'] + '/bin/' + name
        if not os.path.isfile(exe_path):
            raise KeyError
        return exe_path
    except KeyError:
        logger.debug(
            "Unable to find executable in $AMBERHOME/bin. Looking in PATH.")

    exe_path = distutils.spawn.find_executable(name)
    if exe_path is None:
        msg = "Unable to find %s executable, please make sure this is" \
            " present in your PATH or $AMBERHOME/bin." % name
        logger.error(msg)
        raise simulationobjects.SetupError(msg)

    return exe_path


def run_antechamber(lig, charge, resnam=None):
    """
    Wrapper for antechamber with default parameters and AM1-BCC charges

    Parameters
    ----------
    lig : PDBFile or string
      the ligand to run Antechamber on
    charge : int
      the net charge of the ligand
    resnam : string, optional
      the residue name of the ligand

    Returns
    -------
    string
      the filename of the created prepi file
    """

    logger.debug("Running run_antechamber with arguments: ")
    logger.debug("\tlig    = %s" % lig)
    logger.debug("\tcharge = %d" % charge)
    logger.debug("\tresnam = %s" % resnam)
    logger.debug(
        "This will generate an Amber prepi file with AM1-BCC and GAFF atom types"
    )

    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig, six.string_types):
        name = lig
    else:
        name = lig.name

    if resnam is None:
        resnamstr = ""
    else:
        resnamstr = "-rn " + resnam

    ante_exe = _get_executable_path('antechamber')

    # Remove the extension from the filename
    out_name = os.path.splitext(name)[0]
    cmd = '%s -i %s -fi pdb -o %s.prepi -fo prepi -c bcc -nc %d %s -pf y' % (
        ante_exe, name, out_name, charge, resnamstr)
    _run_program('antechamber', cmd)
    subprocess.call("rm sqm.in sqm.out sqm.pdb", shell=True)
    return "%s.prepi" % out_name


def run_parmchk(lig):
    """
    Wrapper for parmcheck with default parameters

    Parameters
    ----------
    lig : PDBFile or string
      the ligand to run Parmcheck on

    Returns
    -------
    string
      the filename of the created frcmod file
    """

    logger.debug("Running run_parmchk with arguments: ")
    logger.debug("\tlig = %s" % lig)
    logger.debug(
        "This will generate an Amber frcmod file with additional parameters")

    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig, six.string_types):
        name = lig
    else:
        name = lig.name

    parm_exe = _get_executable_path('parmchk')

    # Remove the extension from the filename
    out_name = os.path.splitext(name)[0]
    cmd = '%s -i %s.prepi -f prepi -o %s.frcmod' % (parm_exe, out_name,
                                                    out_name)
    _run_program('parmchk', cmd)
    return "%s.frcmod" % out_name


def get_arg_parser():

    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to run antechamber and parmchk"
                    " for a series of PDB-files")
    parser.add_argument(
        '-f', '--files', nargs="+", help="the name of the PDB-files")
    parser.add_argument(
        '-n', '--name', help="the name of the solute", default="UNK")
    parser.add_argument(
        '-c',
        '--charge',
        nargs="+",
        type=float,
        help="the net charge of each PDB-file")
    return parser


if __name__ == "__main__":
    import simulationobjects
    args = get_arg_parser().parse_args()
    # Setup the logger
    logger = simulationobjects.setup_logger("ambertools_py.log")

    for i, filename in enumerate(args.files):
        if args.charge is None or i >= len(args.charge):
            charge = 0.0
        else:
            charge = args.charge[i]
        run_antechamber(filename, charge, args.name)
        run_parmchk(filename + "P")
