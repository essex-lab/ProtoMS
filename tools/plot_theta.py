# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Program to plot the theta distribution, resulting from a JAWS
stage one simulation

Can be executed from the command line as a stand-alone program
"""

import logging
import matplotlib
import numpy as np
import os
from protomslib import simulationobjects

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pylab as plt

logger = logging.getLogger('protoms')


def solname_to_id(molecule, restart='restart'):
    """
  Search for the residue name of a solute in the
  restart file and find out the correspondent
  solute id (or ids, if there are several copies)

  Parameters
  ----------
  molecule : string
    the residue name of the JAWS
    molecule
  restart : string, optional
    the name of the restart file

  Returns
  -------
  list :
    the list of solute / gcsolute ids that
    correspond to that residue name
  string :
    the molecule type. Either "solute", "gcsolute"
    or "both"
  """

    restart_obj = simulationobjects.RestartFile(filename=restart)
    sol_list = [
        key for key in restart_obj.solutesdic
        if restart_obj.solutesdic[key].lower() in molecule.lower()
    ]
    gcsol_list = [
        key for key in restart_obj.gcsolutesdic
        if restart_obj.gcsolutesdic[key].lower() in molecule.lower()
    ]

    if not sol_list:
        return gcsol_list, "gcsolute"
    elif not gcsol_list:
        return sol_list, "solute"
    else:
        return [gcsol_list, sol_list], "both"


def find_solute_theta(molecule, restart='restart', results='results', skip=0):
    """
  Given a solute / gcsolute residue name, find
  and store its corresponding values of theta through
  the simulation.

  Parameters
  ----------
  molecule : string
    the residue name of the JAWS
    molecule
  restart : string, optional
    the name of the restart file
  results : string, optional
    the name of the results file
  skip : int, optional
    the number of snapshots to skip

  Returns
  -------
  dict :
    the theta values per solute id
    in snapshot order
  """

    id_list, mol_type = solname_to_id(molecule, restart=restart)
    results_obj = simulationobjects.ResultsFile()
    results_obj.read(filename=results, skip=skip)
    thetas = dict([[i, []] for i in id_list])
    for snap in results_obj.snapshots:
        if mol_type is "solute":
            res_thetavals = snap.thetasolvals
        elif mol_type is "gcsolute":
            res_thetavals = snap.thetavals
        for ind, res_theta in enumerate(res_thetavals):
            if ind + 1 in thetas.keys():
                thetas[ind + 1].append(res_theta)
    return thetas


def plot_dist(dict_th, plotname="theta_dist", xlab="theta"):
    """
  Plot a distribution of dictonary values
  in form of lists

  Parameters
  ----------
  dict_th : dict
    a dictionary with several keys
    which 'values' are lists of floats
    to plot in a histogram
  plotname : string, optional
    the beginning of the file name
    for the plots produced
  xlab : string, optional
    what to write in the x axes label
    of the histogram

  Return
  ------
  Nothing
    the plots will be saved in png files
  """

    if xlab == "theta":
        xlab = r"$\theta$"

    def set_axes():
        plt.xlabel(xlab, fontsize=20)
        plt.ylabel("Frequency", fontsize=15)

    thetas_values = np.array(dict_th.values(), float)

    n, bins, patches = plt.hist(
        np.reshape(thetas_values, newshape=-1),
        bins=100,
        facecolor=simulationobjects.color(0),
        histtype='stepfilled')
    set_axes()
    plt.savefig(plotname + "_all.png")
    plt.clf()

    n, bins, patches = plt.hist(
        np.transpose(thetas_values), bins=100, histtype="step")
    set_axes()
    plt.savefig(plotname + "_each.png")


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to plot the theta distribution of a given"
                    " molecule, result from a JAWS simulation")
    parser.add_argument(
        '-r',
        '--results',
        help="the name of the results file. Deafult='results'",
        default='results')
    parser.add_argument(
        '-s',
        '--restart',
        help="the replica values to plot. Default='restart'",
        default='restart')
    parser.add_argument(
        '-m',
        '--molecule',
        help="the residue name of the JAWS molecule. Default='WAT'",
        default="WAT")
    parser.add_argument(
        '-p',
        '--plotname',
        help="the start of the filename for the plots generated. "
             "Default='theta_dist'",
        default="theta_dist")
    parser.add_argument(
        '--skip',
        help="the number of results snapshots to skip, Default = 0",
        default="0")
    return parser


#
# If this is run from the command-line
#
if __name__ == '__main__':

    args = get_arg_parser().parse_args()

    try:
        skip_steps = int(args.skip)
    except (ValueError, TypeError):
        simulationobjects.SetupError(
            "Snapshots to skip needs to be an integer. The argument %s"
            " could not be interpreted." % args.skip)

    thetas_dic = find_solute_theta(args.molecule, args.restart, args.results,
                                   skip_steps)

    plot_dist(thetas_dic, plotname=args.plotname)
