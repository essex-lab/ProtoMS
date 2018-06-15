# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Program to analyze and plot data series

This module defines the following public functions:
find_equilibration
stats_inefficiency
maximize_samples
running
moving
parse_series
parse_files
plot_series
write_series

Can be executed from the command line as a stand-alone program
"""

from __future__ import print_function
import logging
import matplotlib
import numpy as np
import os
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from protomslib import simulationobjects


if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pylab as plt


logger = logging.getLogger('protoms')

EQUIL_LIMIT = 10

################################
# Data series analysis routines
################################


def find_equilibration(x, y, atleast=EQUIL_LIMIT, threshold=0.05, nperm=0):
    """
  Find the equilibration time of a data series

  Parameters
  ----------
  x : Numpy array
    the x data
  y : Numpy array
    the y data
  atleast : int, optional
    this is the smallest number of snapshots considered to be production
  threshold : float, optional
    the confidence level
  nperm : int, optional
    if greater than zero, a permutation test is peformed
    with nperm synthetic permutations

  Returns
  -------
  int
    the length of the equilibration period
  """
    x, y = np.array(x), np.array(y)
    for i in range(0, y.shape[0] - atleast):
        if nperm == 0:  # Performs an assymptotic test, fast
            tau, p = kendalltau(x[i:], y[i:])
            if p > threshold:
                return i
        else:  # Performs a rigorous permutation test, slow
            rhos = np.zeros(nperm)
            rho0, p = spearmanr(x[i:], y[i:])
            for j in range(nperm):
                rhos[j], p = spearmanr(x[i:], np.random.permutation(y[i:]))
            ncnt = (rho0 * rho0 < rhos * rhos).sum()
            if float(ncnt) / float(nperm) > threshold:
                return i
    return y.shape[0] - atleast


def stat_inefficiency(y):
    """
  Calculates the statistical inefficiency, g.

  g = 1 + 2*tau, where tau is the autocorrelation time

  Parameters
  ----------
  y : Numpy array
    the data

  Returns
  -------
  float
    g, the inefficiency
    or None if variance is too small
  float
    neff, the number of samples 
    or None if variance is too small
  """
    n = len(y)
    dy = y - y.mean()
    vary = (dy * dy).mean()
    if vary < 0.1E-6:
        return None, None
    tmax = int(np.round(n * 0.9))
    g = 1.0
    for t in range(1, tmax + 1):
        c = np.sum(dy[0:n - t] * dy[t:n]) / (float(n - t) * vary)
        if c < 0.0 and t > 10:
            break
        g = g + 2.0 * c * (1.0 - float(t) / float(n))
    if g < 1.0:
        g = 1.0
    neff = float(n) / float(g)
    return g, neff


def maximize_samples(y, atleast=EQUIL_LIMIT):
    """
  Find the minimum of the statistical inefficiency by varying
  the equlibration period, i.e. maximizing the number of
  uncorrelated samples

  Parameters
  ----------
  y : Numpy array
    the data
  atleast : int, optional
    this is the smallest number of snapshots considered to be production

  Returns
  -------
  int
    the snapshot at which g is maximum
  float
    maximum of g
  int
    the number of uncorrelated samples
  """
    n = y.shape[0]
    gt = np.zeros([n - atleast])
    neff = np.zeros([n - atleast])
    for i in range(0, n - atleast):
        gt[i], neff[i] = stat_inefficiency(y[i:])
        # If stat_inefficiency return None, we set g and neff to 0
        if np.isnan(gt[i]):
            gt[i], neff[i] = 0.0, 0.0

    neff_max = neff.max()
    t = neff.argmax()
    return t, gt[t], neff_max


def running(data):
    """
  Returns a running average of the data

  Value at position i of the returned data is the average of the
   input data from point zero to point i

  Parameters
  ----------
  data : NumpyArray
    the input data, should be 1D

  Returns
  -------
  NumpyArray
    the averaged data
  """
    running = np.zeros(data.shape)
    for i in range(data.shape[0]):
        running[i] = data[:i + 1].mean()
    return running


def moving(data, window):
    """
  Returns a moving/rolling average of the data

  Parameters
  ----------
  data : NumpyArray
    the input data, should be 1D
  window : int
    the window size

  Returns
  -------
  NumpyArray
    the averaged data
  """
    weights = np.repeat(1.0, window) / float(window)
    return np.convolve(data, weights, 'valid')


###################
# Parsing routines
###################


def parse_series(series, results):
    """
  Extract data series and label from a results file based on a label

  Parameters
  ----------
  series : list of string
    the series to extract
  results : SnapshotResults object
    all the results

  Returns
  -------
  list of NumpyArrays
    the data series
  list of string
    the labels for the series
  """

    def parse_compound_series(series, results):
        """
        Parse a compound series, i.e. internal or interaction energies
        """
        cols = series.strip().split("/")
        s2a = {"intra": "internal_energies", "inter": "interaction_energies"}
        attr = s2a[cols[0]]
        elabel = cols[1]
        etype = cols[2]

        label = elabel + "/" + etype + " [kcal/mol]"
        if etype[-1] == "b":
            eattr = "back"
            etype = etype[:-1]
        elif etype[-1] == "f":
            eattr = "forw"
            etype = etype[:-1]
        else:
            eattr = "curr"

        if hasattr(results, attr):
            dict = getattr(results, attr)
            if elabel in dict:
                for e in dict[elabel]:
                    if e.type.lower() == etype:
                        return getattr(e, eattr), label
        return None

    def parse_energyresults(series, results):
        """
    Parse a EnergyResults object
    """
        if series[-1] == "b":
            eattr = "back"
            etype = series[:-1]
        elif series[-1] == "f":
            eattr = "forw"
            etype = series[:-1]
        else:
            eattr = "curr"
            etype = series
        label = series + " [kcal/mol]"
        if hasattr(results, etype):
            return getattr(getattr(results, etype), eattr), label
        else:
            return None

    def parse_feenergy(series, results):
        """
    Parse the feenergy attribute
    """
        if not hasattr(results, "feenergies"): return None

        cols = series.lower().strip().split("/")
        if len(cols) == 1:
            ys = [
                results.feenergies[l]
                for l in sorted(results.feenergies.keys())
            ]
            labels = [
                "energy at %.3f [kcal/mol]" % l
                for l in sorted(results.feenergies.keys())
            ]
            return ys, labels
        else:
            ys = [
                results.feenergies[float(l)] for l in cols[1].split(",")
                if float(l) in results.feenergies
            ]
            labels = [
                "energy at %.3f [kcal/mol]" % float(l)
                for l in cols[1].split(",") if float(l) in results.feenergies
            ]
            if len(ys) == 0:
                return None
            else:
                return ys, labels

    # Setup parser for special keywords
    special_parse = {}
    for attr in ["total", "capenergy", "extraenergy"]:
        special_parse[attr] = parse_energyresults
        special_parse[attr + "f"] = parse_energyresults
        special_parse[attr + "b"] = parse_energyresults
    special_parse["intra"] = parse_compound_series
    special_parse["inter"] = parse_compound_series
    special_parse["feenergy"] = parse_feenergy

    # Setup units for special keywords, if not defined here,
    # the unit will be kcal/mol
    special_units = {}
    special_units["solventson"] = ""
    special_units["lambdareplica"] = ""
    special_units["globalreplica"] = ""

    ys = []
    labels = []
    for s in args.series:
        resp = None
        attr = s.strip().split("/")[0]
        if attr in special_parse:
            resp = special_parse[attr](s, results)
        else:
            if hasattr(results, s):
                if s in special_units:
                    l = s + special_units[s]
                else:
                    l = s + " [kcal/mol]"
                resp = getattr(results, s), l
        if resp is not None:
            y, label = resp
            if isinstance(y, list):
                ys.extend(y)
                labels.extend(label)
            else:
                ys.append(y)
                labels.append(label)
        else:
            print("Skipping non-exisiting series %s" % s)
    return ys, labels


def parse_files(filenames, series, all_results, ys, labels, **kwargs):
    """
  Parse multiple results file

  If a data serie in series can be calculated over multiple
  input files, the data series is replaced by some calculation
  over all input files

  Parameters
  ----------
  filenames : list of strings
    the files to parse
  series : list of strings
    the series selected
  all_results : list of SnapshotResults
    all the results loaded from all input files, needs to
     be initialized before calling this routine
  ys : list of NumpyArrays
    all the series
  labels : list of strings
    labels for all the data series
  **kwargs : dictionary
    additional parameters, passed directly to the different parsers

  Returns
  -------
  bool
    whether the series has been modified, i.e. if any parser was found
  """

    def parse_multi_gradients(results, attr, **kwargs):
        """
    Routine to return parse multiple gradients series,
    and computing running average of dG using trapezium
    """
        lam = [r.lam[0] for r in results]
        lam = np.array(lam)

        gradients = [getattr(r, attr) for r in results]
        gradients = np.array(gradients)

        dg = np.zeros(gradients.shape[1])

        # Setup moving averages
        domoving = "moving" in kwargs and kwargs["moving"] is not None
        if domoving:
            d = int(np.ceil(kwargs["moving"] / 2.0))
        else:
            d = 0

        for i in range(d, gradients.shape[1] - d):
            if domoving:
                dg[i] = np.trapz(
                    gradients[:, i - d:i + d + 1].mean(axis=1), lam)
            else:
                dg[i] = np.trapz(gradients[:, i:].mean(axis=1), lam)

        # Remove start and begininng of series
        if domoving:
            dg = dg[d:dg.shape[0] - d]

        return dg, "dG [kcal/mol]"

    # Main routine

    # Initialize arrays
    if len(ys) == 0:
        ys = [None] * len(series)
        labels = [None] * len(series)

    # Initialize parsing routines
    parse_multi = {}
    parse_multi["gradient"] = parse_multi_gradients
    parse_multi["agradient"] = parse_multi_gradients

    # Check that some parsing routines for the selected series exists
    found = False
    for i, s in enumerate(series):
        attr = s.strip().split("/")[0]
        if attr in parse_multi:
            found = True
            break

    if found:
        # Open all results files
        for f in filenames:
            results_file = simulationobjects.ResultsFile()
            results_file.read(filename=f)
            all_results.append(results_file.make_series())
        # Replace the series with a series created from all input files
        for i, s in enumerate(series):
            attr = s.strip().split("/")[0]
            if attr in parse_multi:
                ys[i], labels[i] = parse_multi[attr](all_results, s, **kwargs)
    return found


#################################################
# Routines to plot and write data series to disc
#################################################


def _label0(label):
    """
  Removes the unit from a label
  """
    if label.find("[") == -1:
        return label

    l = "_".join(label.split()[:-1])
    l = l.replace("/", "_")
    return l


def plot_series(ys, yprop, labels, offset, plotkind, outprefix):
    """
  Plot a series

  Parameters
  ----------
  ys : list of Numpy arrays
    the data series
  yprop : list of dictionary
    statistical properties of the data
  labels : list of strings
    the labels for the data series
  offset : int
    the offset of x from 1
  plotkind : string
    the type of plot to create, can be either sep, sub,
    single, single_first0 or single_last0
  outprefix : string
    the prefix of the created png-files
  """
    if plotkind not in ["sep", "sub", "single",
                        "single_first0", "single_last0"]:
        plotkind == "sep"
    if len(ys) == 1 or plotkind is None:
        plotkind = "single"  # For a single series, there is only one kind
    ncols = int(np.ceil(len(ys) / 2.0))

    ys = np.array(ys)
    x = np.arange(1, ys.shape[1] + 1) + offset

    # Subtract either first or last point from each series
    if plotkind.startswith("single_"):
        for i in range(ys.shape[0]):
            if plotkind == "single_first0":
                ys[i, :] = ys[i, :] - ys[i, 0]
            else:
                ys[i, :] = ys[i, :] - ys[i, -1]

    if plotkind != "sep":
        currfig = plt.figure(1)

    if plotkind not in ["sep", "sub"]:
        ymax = ys.max() + 0.1 * (ys.max() - ys.min())
        ymin = ys.min() - 0.1 * (ys.max() - ys.min())

    for i, (y, prop, label) in enumerate(zip(ys, yprop, labels)):
        if plotkind == "sep":
            currfig = plt.figure(i + 1)
        if plotkind == "sub":
            ax = currfig.add_subplot(2, ncols, i + 1)
        else:
            ax = currfig.gca()

        if plotkind in ["sep", "sub"]:
            ymax = y.max() + 0.1 * (y.max() - y.min())
            ymin = y.min() - 0.1 * (y.max() - y.min())

        ax.plot(x, y, label=label, color=simulationobjects.color(i))
        ax.plot(
            [prop["equil"], prop["equil"]], [ymin, ymax],
            '--',
            color=simulationobjects.color(i))

        if not plotkind.startswith("single"):
            ax.set_ylim([ymin, ymax])
            ax.set_xlabel("Snapshot")
            ax.set_ylabel(label)
        if plotkind == "sep":
            currfig.savefig(
                outprefix + "_" + _label0(label) + ".png", format="png")

    if plotkind.startswith("single"):
        ax.set_ylim([ymin, ymax])
        ax.set_xlabel("Snapshot")
        if ys.shape[0] == 1:
            ax.set_ylabel(labels[0])
        else:
            ax.set_ylabel("Data")
            ax.legend()

    if plotkind != "sep":
        currfig.savefig(outprefix + ".png", format="png")
        if "DISPLAY" in os.environ and os.environ["DISPLAY"] != "":
            currfig.show()
            print("\nType enter to quit\n>", end="")
            raw_input()


def write_series(ys, yprop, labels, offset, filekind, outprefix):
    """
  Write a series to disc

  Parameters
  ----------
  ys : list of Numpy arrays
    the data series
  yprop : list of dictionary
    statistical properties of the data
  labels : list of strings
    the labels for the data series
  offset : int
    the offset of x from 1
  filekind : string
    the type of file to write, can be either sep or single
  outprefix : string
    the prefix of the created files
  """
    if len(ys) == 1 or filekind is None:
        filekind = "sep"  # For a single series, there is only one kind
    if filekind not in ["sep", "single"]:
        # For files the kinds sub, single, single_first0 and
        # single_last0 all means the same
        filekind = "single"

    if filekind == "sep":
        for y, prop, label in zip(ys, yprop, labels):
            with open(outprefix + "_" + _label0(label) + ".dat", "w") as f:
                f.write("#Data for %s\n" % label)
                f.write("#Equilibration time: %d\n" % prop["equil"])
                f.write("#Production contain %d data (g=%.3f)\n" %
                        (prop["neff"], prop["g"]))
                f.write("#Number of samples are maximized at %d with g=%.3f"
                        " and samples=%d\n" % (
                            prop["t_opt"], prop["g_min"], prop["neff_max"]))
                for i, yi in enumerate(y):
                    f.write("%d %.5f\n" % (i + 1 + offset, yi))
    else:
        with open(outprefix + ".dat", "w") as f:
            f.write("#Data for %s\n" % "\t".join(labels))
            f.write(
                "#Equilibration time: %s\n" % "\t".join("%d" % prop["equil"]
                                                        for prop in yprop))
            f.write(
                "#Independent samples: %s\n" % "\t".join("%d" % prop["neff"]
                                                         for prop in yprop))
            f.write(
                "#Production period g: %s\n" % "\t".join("%d" % prop["g"]
                                                         for prop in yprop))
            f.write("#Snapshot when independent samples are optimized: %s\n" %
                    "\t".join("%d" % prop["t_opt"] for prop in yprop))
            f.write("#Minimum g: %s\n" % "\t".join("%d" % prop["g_min"]
                                                   for prop in yprop))
            f.write("#Maximum number of independent samples: %s\n" % "\t".join(
                "%d" % prop["neff_max"] for prop in yprop))
            for i in range(ys[0].shape[0]):
                f.write("%d %s\n" % (i + 1 + offset, "\t".join("%.5f" % y[i]
                                                               for y in ys)))


#################################
# Wizard routines for user input
#################################


def _select_series(results):
    """
  Prompts the user for series to plot and write

  Parameters
  ----------
  results : SnapshotResults object
    all the results

  Returns
  -------
  list of string
    the selected series
  """
    selection = []

    print("Select one or several of the following series:")
    print("----------------------------------------------")
    singles = []
    for attr in ["backfe", "forwfe", "gradient", "agradient",
                 "lambdareplica", "solventson", "globalreplica"]:
        if hasattr(results, attr):
            singles.append(attr)
    for attr in ["total", "capenergy", "extraenergy"]:
        if hasattr(results, attr):
            singles.append(attr + "[b|f]")
    if len(singles) > 0:
        print("Single valued series: ")
        for i in range(0, len(singles), 5):
            print("".join("%-15s" % s for s in singles[i:i + 5]))
        print("(type e.g. total)")
    if hasattr(results, "internal_energies"):
        print("\nInternal energies:")
        for elabel in results.internal_energies:
            print("%-10s : %s" % (elabel, ", ".join(
                e.type.lower() for e in results.internal_energies[elabel])))
        print("(type e.g. intra/protein1/sum)")
    if hasattr(results, "interaction_energies"):
        print("\nInteraction energies:")
        for elabel in results.interaction_energies:
            print("%-20s : %s" % (elabel, ", ".join(
                e.type.lower() for e in results.interaction_energies[elabel])))
        print("(type e.g. inter/solvent-solvent/sum)")
    if hasattr(results, "feenergies"):
        print("\nEnergies at lambda-values: %s" % (", ".join(
            "%.3f" % l for l in sorted(results.feenergies.keys()))))
        print("(type e.g. feenergy or feenergy/0.000,1.000)")

    print("\n> ", end="")
    instr = raw_input()
    while len(instr) > 0:
        selection.append(instr)
        print("> ", end="")
        instr = raw_input()
    return selection


def _select_plot():
    """
  Prompt the user to select a plot kind

  Returns
  -------
  string
    the plot/write kind
  """
    print("\nHow do you want to plot/write the multiple series?")
    print("1) Separate plots (default)")
    print("2) Sub-plots")
    print("3) Single plot")
    print("4) Single plot + subtract first snapshot")
    print("5) Single plot + subtract last snapshot")
    print("\n> ", end="")
    instr = raw_input()
    if instr == "2":
        return "sub"
    elif instr == "3":
        return "single"
    elif instr == "4":
        return "single_first0"
    elif instr == "5":
        return "single_last0"
    else:
        return "sep"


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to analyze and plot a time series")
    parser.add_argument(
        '-f',
        '--file',
        nargs="+",
        help="the name of the file to analyse. Default is results. ",
        default=["results"])
    parser.add_argument(
        '-o',
        '--out',
        help="the prefix of the output figure. Default is series. ",
        default="series")
    parser.add_argument(
        '-s', '--series', nargs="+", help="the series to analyze")
    parser.add_argument(
        '-p',
        '--plot',
        choices=["sep", "sub", "single", "single_first0", "single_last0"],
        help="the type of plot to generate for several series")
    parser.add_argument(
        '--nperm',
        type=int,
        help="if larger than zero, perform a permutation test to"
             " determine equilibration, default=0",
        default=0)
    parser.add_argument(
        '--threshold',
        type=float,
        help="the significant level of the equilibration test, default=0.05",
        default=0.05)
    parser.add_argument(
        '--average',
        action='store_true',
        help="turns on use of running averaging of series",
        default=False)
    parser.add_argument(
        '--moving',
        type=int,
        help="turns on use of moving averaging of series, default=None")
    return parser


#
# If this is run from the command-line
#
if __name__ == '__main__':

    args = get_arg_parser().parse_args()

    # Mutually exclusive options
    if args.average:
        args.moving = None

    if args.moving is not None:
        offset = int(np.ceil(args.moving / 2.0))
    else:
        offset = 0

    # Read the input file from disc
    if len(args.file) == 2 and os.path.isdir(args.file[0]):
        restem = args.file.pop(1)
        filename = [args.file[0], restem]
    else:
        filename = args.file[0]
    results_file = simulationobjects.ResultsFile()
    results_file.read(filename=filename)
    results = results_file.make_series(
    )  # This puts all data into Numpy arrays

    # Select which series to plot
    if args.series is None:
        args.series = _select_series(results)
    if len(args.series) == 0:
        print("No series selected! Exiting")
        quit()

    # Parse the user selected series into NumpyArray with the data
    # and an appropriate label
    ys, labels = parse_series(args.series, results)
    if len(ys) == 0:
        print("No series was parsed from the selection! Exiting")
        quit()

    # If user has given several results files
    parsed_files = False
    if len(args.file) > 1:
        all_results = [results]
        parsed_files = parse_files(
            args.file[1:],
            args.series,
            all_results,
            ys,
            labels,
            moving=args.moving)

    # Calculates running/moving averages for series
    if not parsed_files and args.average:
        for i, y in enumerate(ys):
            ys[i] = running(y)
    elif not parsed_files and args.moving is not None:
        for i, y in enumerate(ys):
            ys[i] = moving(y, args.moving)

    # Compute the equilibration time and other properties for
    # each series to plot
    yprop = [{
        "equil": 0,
        "g": 0,
        "neff": 0,
        "t_opt": 0,
        "neff_max": 0,
        "g_min": 0
    } for y in ys]
    x = np.arange(1, ys[0].shape[0] + 1) + offset
    print()
    for i, (y, prop, label) in enumerate(zip(ys, yprop, labels)):
        if len(y) <= EQUIL_LIMIT:
            print("%i snapshots or less found in %s, will not evaluate "
                  "equilibration" % (EQUIL_LIMIT, _label0(label)))
            continue
        prop["equil"] = find_equilibration(
            y, x, nperm=args.nperm, threshold=args.threshold) + offset
        if prop["equil"] == EQUIL_LIMIT:
            print("No point of equilibration found for %s" % _label0(label))
            prop["equil"] = 0
            continue

        print("Equilibration found at snapshot %d for %s, value=%.3f" %
              (prop["equil"], _label0(label), ys[i][prop["equil"]]))
        prop["g"], prop["neff"] = stat_inefficiency(y[prop["equil"]:])
        if prop["g"] is not None:
            print("\tThis production part is estimated to contain %d "
                  "uncorrelated samples (g=%.3f)" % (prop["neff"], prop["g"]))
            prop["t_opt"], prop["g_min"], prop["neff_max"] = \
                maximize_samples(y)
            prop["t_opt"] = prop["t_opt"] + offset
            print("\tThe number of samples is maximized at %d, g=%.3f and the "
                  "number of uncorrelated samples is %d" % (
                      prop["t_opt"], prop["g_min"], prop["neff_max"]))

    # Select what kind of plot to make for multiple series
    if len(ys) > 1 and args.plot is None:
        args.plot = _select_plot()

    # Plot the series
    plot_series(ys, yprop, labels, offset, args.plot, args.out)

    # Write the series to disc
    write_series(ys, yprop, labels, offset, args.plot, args.out)
