"""Collection of classes to form the basis of a replacement for the current
free energy calculation framework. """

from copy import copy
from glob import glob
from operator import add
import os
import matplotlib
import numpy as np
import pickle
import pymbar
from scipy.integrate import trapz
import warnings
import simulationobjects as sim
from free_energy_argument_parser import FEArgumentParser

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['lines.linewidth'] = 3


class Series(object):
    """Object for storing data series, i.e. some value that changes
    as a function of some coordinate. Designed to act as a base class
    for PMF classes.
    """
    def __init__(self, coordinate, *args):
        """Parameters
        ----------
        coordinate : sequence of numbers
            The values of the coordinate for the PMF
        *args : An arbitrary number of iterables containing numbers
                or Quantity objects
            All arguments after coordinate are iterated simultaneously,
            the first entry of each is used to construct the first Quantity
            object in self.values and so on."""
        self.coordinate = coordinate
        self.values = []
        for dat in zip(*args):
            self.values.append(Quantity.fromData(dat))

    def __add__(self, other):
        if list(self.coordinate) != list(other.coordinate):
            raise ValueError(
                'Cannot add PMF instances with different coordinate values')
        new_pmf = copy(self)
        new_pmf.values = [s + o for s, o in zip(self.values, other.values)]
        return new_pmf

    def __neg__(self):
        other = copy(self)
        other.values = [-val for val in other.values]
        return other

    def __iter__(self):
        return iter(self.values)

    def plot(self, axes, show_error=True, fmt='-', xlabel='Lambda Value',
             ylabel='Free Energy (kcal)', **kwargs):
        """Plot this PMF onto the provided figure axes.

        Parameters
        ----------
        axes : matplotlib Axes object
            axes onto which the figure will be drawn
        xlabel : string
            label for x-axis
        ylabel : string
            label for y-axis
        **kwargs :
            additional keyword arguments to be passed to axes.plot
        """
        y = np.array([fe.value for fe in self.values])
        err = np.array([fe.error for fe in self.values])

        line = axes.plot(self.coordinate, y, fmt, **kwargs)[0]
        if show_error:
            axes.plot(self.coordinate, y+err, '--',
                      linewidth=1, color=line.get_color())
            axes.plot(self.coordinate, y-err, '--',
                      linewidth=1, color=line.get_color())
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        return line


class PMF(Series):
    """A Potential of Mean Force, describing a free energy profile
    as a function of some coordinate.
    """

    @property
    def dG(self):
        """The free energy difference at the PMF end points"""
        return self.values[-1] - self.values[0]


class BaseResult(object):
    """A base class for Result objects."""
    def __init__(self, *args, **kwargs):
        """Parameters
        ----------
        *args : an arbitrary number of lists of PMF objects
            Each list should contain the PMF objects corresponding to
            multiple repeats within a calculation leg. The free energy
            estimates from multiple lists will be summed together."""
        self.data = args
        if len({tuple(pmf.coordinate)
                for dat in self.data for pmf in dat}) > 1:
            raise ValueError("All data must use the same lambda values")
        # try:
        #     self.pmf_class = kwargs['pmf_class']
        # except KeyError:
        #     self.pmf_class = PMF
        self.pmf_class = PMF

    def __add__(self, other):
        return self.__class__(*(self.data + other.data),
                              pmf_class=self.pmf_class)

    def __sub__(self, other):
        return self + -other

    def __neg__(self):
        return self.__class__(*[[-pmf for pmf in dat] for dat in self.data],
                              pmf_class=self.pmf_class)


class Result(BaseResult):
    """The result of a free energy calculation. Contains multiple PMF
    objects that are combined together to give a meaningful free energy
    difference.
    """
    @property
    def lambdas(self):
        """Lambda values at which the calculation has been performed"""
        try:
            return self.data[0][0].coordinate
        except IndexError:
            raise ValueError("This result does not contain any data.")

    @property
    def dG(self):
        """The free energy difference at the PMF end points"""
        data_dGs = [[pmf.dG for pmf in dat] for dat in self.data]
        FEs = [Quantity.fromData(dG) for dG in data_dGs]
        return sum(FEs, Quantity(0., 0.))

    @property
    def pmf(self):
        """The potential of mean force for this result"""
        pmfs = [self.pmf_class(self.lambdas, *dat) for dat in self.data]
        return reduce(add, pmfs)


class Quantity(object):
    """Stores both the value and error associated with a
    free energy estimate.
    """
    def __init__(self, value, error):
        """
        Parameters
        ----------
        value: float
            the magnitude of the free energy difference
        error: float
            one standard error associated with value
        """
        self.value = value
        self.error = error

    @staticmethod
    def fromData(data):
        """Alternative initialiser that uses of a set of data
        points rather than providing an explicit value and error.
        Returns a Quantity object with a value and error determined
        from the mean of and standard error of the provided data. If
        Quantity objects are provided as data their individual
        error attributes are ignored

        Parameters
        ----------
        data: sequence of numbers or Quantity objects
           data from which to determine value and error
        """

        # suppress warnings that arise when data is an empty list
        # NaN is the output and this is desired
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            try:
                return Quantity(np.mean(data), np.std(data)/len(data)**0.5)
            except TypeError:
                values = [dat.value for dat in data]
                return Quantity(np.mean(values),
                                np.std(values)/len(values)**0.5)

    def __add__(self, other):
        return Quantity(self.value + other.value,
                        (self.error**2 + other.error**2)**0.5)

    def __sub__(self, other):
        return Quantity(self.value - other.value,
                        (self.error**2 + other.error**2)**0.5)

    def __neg__(self):
        return Quantity(-self.value, self.error)

    def __eq__(self, other):
        return (self.value == other.value) and (self.error == other.error)

    def __str__(self):
        return "%9.4f +/- %.4f" % (self.value, self.error)

    def __repr__(self):
        return "<Quantity: value=%.4f error=%.4f>" % (self.value, self.error)


class Estimator(object):
    """Base class for free energy estimators."""
    results_name = 'results'
    # subdir_glob = './'
    pmf_class = PMF
    result_class = Result

    def __init__(self, lambdas, subdir_glob="./", **kwargs):
        """Parameters
        ----------
        lambdas: list of numbers
            lambda values used in a calculation
        """
        self.data = []
        self.lambdas = lambdas
        self.subdir_glob = subdir_glob

    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular lambda value
        """
        pass

    def calculate(self, temp=300.):
        """Calculate the free energy difference and return a PMF object.

        Parameters
        ----------
        temp: float, optional
              temperature of calculation
        """
        pass

    def __getitem__(self, val):
        """Return a new class instance with series[val] applied to each
        individual data series.
        """
        new_est = copy(self)
        new_est.data = []
        # add data series to the new estimator that have been sliced by val
        # want to always apply slice to last dimension, so transpose array
        # apply slice to first dimension and then transpose back
        for dat in self.data:
            reordered_dat = dat.T[val]
            new_est.data.append(reordered_dat.T)
        return new_est

    def __len__(self):
        """Return the number of data points contained in the data series."""
        lengths = [dat.shape[-1] for dat in self.data]
        if len(lengths) == 0:
            return 0
        if len(set(lengths)) > 1:
            raise Exception("Found data entries of different lengths.")
        return lengths[0]

    def subset(self, low_bound=0., high_bound=1., step=1):
        """Return a new class instance containing truncated data series.

        Parameters
        ----------
        low_bound: float, default=0.
            Starting position of truncated data series.
        high_bound: float, default=1.
            Finishing position of truncated data series.
        step: int, default=1
            Subsampling frequency of data. 1=all data points,
            2=every other data point, etc.
        """
        if not(0. <= low_bound < high_bound <= 1.):
            raise ValueError("Both bounds must be in the range [0,1] and "
                             "low_bound must be less than high_bound.")

        low_ind = int(len(self)*low_bound)
        high_ind = int(len(self)*high_bound)
        if low_ind == high_ind:
            raise ValueError("Specified bounds will return estimator "
                             "containing no data.")

        return self[low_ind:high_ind:step]


class TI(Estimator):
    """Estimate free energy differences using Thermodynamic Integration."""
    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular lambda value
        """
        self.data.append((series.forwfe-series.backfe) /
                         (series.lamf[0]-series.lamb[0]))
        return len(self.data[-1])

    def calculate(self, temp=300.):
        """Calculate the free energy difference and return a PMF object.

        Parameters
        ----------
        temp: float, optional
              temperature of calculation
        """
        gradients = [gradient_data.mean() for gradient_data in self.data]
        pmf_values = [trapz(gradients[:i], self.lambdas[:i])
                      for i in xrange(1, len(self.lambdas) + 1)]
        return PMF(self.lambdas, pmf_values)


class BAR(Estimator):
    """Estimate free energy differences using Bennett's Acceptance Ratio."""
    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular lambda value
        """
        lam = series.lam[0]
        lam_ind = self.lambdas.index(lam)
        lamf = self.lambdas[lam_ind+1] if lam != 1.0 else 1.0
        lamb = self.lambdas[lam_ind-1] if lam != 0.0 else 0.0
        self.data.append(
            np.array([series.feenergies[lam] - series.feenergies[lamb],
                      series.feenergies[lam] - series.feenergies[lamf]]))
        return len(self.data[-1][0])

    def calculate(self, temp=300.):
        """Calculate the free energy difference and return a PMF object.

        Parameters
        ----------
        temp: float, optional
              temperature of calculation
        """
        beta = 1./(sim.boltz*temp)
        pmf_values = [0.0]
        for low_lam, high_lam in zip(self.data, self.data[1:]):
            pmf_values.append(pmf_values[-1] +
                              pymbar.BAR(-low_lam[1]*beta,
                                         -high_lam[0]*beta)[0]/beta)
        return PMF(self.lambdas, pmf_values)


class MBAR(TI):
    """Estimate free energy differences using the Multistate Bennett's
    Acceptante Ratio.
    """
    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular lambda value
        """
        self.data.append(np.array([series.feenergies[lam]
                                   for lam in sorted(series.feenergies)]))
        return len(self.data[-1][0])

    def calculate(self, temp=300.):
        """Calculate the free energy difference and return a PMF object.

        Parameters
        ----------
        temp: float, optional
              temperature of calculation
        """
        beta = 1./(sim.boltz*temp)
        mbar = pymbar.MBAR(np.array(self.data)*beta,
                           [len(dat[0]) for dat in self.data])
        FEs = mbar.getFreeEnergyDifferences(compute_uncertainty=False)[0]/beta
        return PMF(self.lambdas,
                   [FEs[0, i] for i in xrange(len(self.data))])


class FreeEnergyCalculation(object):
    """Class for performing a free energy calculation from one or more
    ProtoMS simulation outputs.
    """

    def __init__(self, root_paths, temperature,
                 estimators=[TI, BAR, MBAR], subdir='',
                 extract_data=True, **kwargs):
        """Parameters
        ----------
        root_paths: a list of lists of strings
            Paths to ProtoMS output directories. Each list of strings specifies
            the output directory(s) that make up multiple repeats of a single
             'leg' of the calculation. The result for each leg is a mean over
            constituent repeats. Results from each leg are then summed together
            to give a total free energy difference for the calculation.
        temperature: float
            Simulation temperature in degrees Kelvin.
        estimators: list of Estimator classes
            Estimator classes to use in evaluating free energies.
        subdir: string
            subdirectory within each lambda folder to search for output
        extract_data: bool
            If True extract data from simulation results in prepartion for
            calculation.
        **kwargs:
            Additional keyword arguments are passed to Estimator classes
            during initialisation.
        """
        self.root_paths = root_paths
        self.temperature = temperature
        # self.subdir = subdir
        self.estimator_classes = estimators

        # data hierarchy -> root_paths[leg][repeat]
        self.paths = []
        self.lambdas = []
        for leg in self.root_paths:
            self.paths.append([])
            self.lambdas.append([])
            for root_path in leg:
                paths = glob(self._path_constructor(root_path))
                paths.sort(key=self._get_lambda)
                self.paths[-1].append(paths)
                self.lambdas[-1].append(map(self._get_lambda, paths))

        self.estimators = {
            estimator: [[estimator(l, subdir_glob=subdir, **kwargs)
                         for l in lams]
                        for lams in self.lambdas]
            for estimator in estimators}
        if extract_data:
            self._extract_series(self.paths)
        self.figures = {}
        self.tables = []

    def _path_constructor(self, root_path):
        """Given a root_path (string) construct a full path
        suitable for globbing to find output directories."""
        return os.path.join(root_path, "lam-*") #, self.subdir)

    def _get_lambda(self, path):
        """Given a path (string) extract the contained lambda value"""
        return float(path.split('/')[-1][4:])

    def _extract_series(self, paths):
        """For each entry of paths (list of strings) open the contained
        results file and extract the data series and supply these
        to each estimator instance."""
        for i, leg in enumerate(paths):
            min_len = 10E20
            for j, repeat in enumerate(leg):
                for path in repeat:
                    # get the names of required data files
                    # for each estimator class
                    result_names = {}
                    for cls in self.estimators:
                        est = self.estimators[cls][i][j]
                        file_path = os.path.join(
                            path,
                            est.subdir_glob,
                            est.results_name)
                        result_names[cls] = glob(file_path)
                        if result_names[cls] == []:
                            raise Exception(
                                "Unable to find ProtoMS output files matching "
                                "the file glob:\n%s" % file_path)
                    # get the unique list of file names so that we only
                    # open each one once for efficiency
                    result_names_uniq = set(v for val in result_names.values()
                                            for v in val)
                    for fname in sorted(result_names_uniq):
                        rf = sim.ResultsFile()
                        rf.read(fname)
                        rf.make_series()
                        for est in self.estimators.values():
                            inst = est[i][j]
                            # only add data from a results file if it is in
                            # list for the associated class
                            if fname in result_names[inst.__class__]:
                                data_len = inst.add_data(rf.series)
                                min_len = data_len \
                                    if data_len < min_len else min_len

                # in cases where calculations terminate prematurely there
                # can be slight differences in the length of data series
                # within a repeat. Here we standardise the length for
                # later convenience
                for est in self.estimators.values():
                    est[i][j] = est[i][j][:min_len]

    def calculate(self, subset=(0., 1., 1)):
        """For each estimator return the evaluated potential of mean force.

        Parameters
        ----------
        subset: tuple of two floats and an int
            specify the subset of data to use in the calculation
        """
        results = {}
        for est, legs in self.estimators.items():
            leg_result = est.result_class()
            for i, leg in enumerate(legs):
                leg_result += est.result_class(
                    [rep.subset(*subset).calculate(self.temperature)
                     for rep in leg])
            results[est] = leg_result
        return results

    def run(self, args):
        """Public method for execution of calculation. This is usually the
        desired way to use execute calculations. Handles logic for
        printing of output tables and saving of output pickles and
        figures. Calls plt.show() to display figures.  Runs
        self._body(args) and uses the return value as the calculation
        result.

        Parameters
        ----------
        args
        """

        results = self._body(args)

        if args.pickle is not None:
            with open(args.pickle, 'w') as f:
                pickle.dump(results, f, protocol=2)

        for table in self.tables:
            print "%s\n" % table

        if args.save_figures is not None:
            # flag provided so save figures
            for key in self.figures:
                if args.save_figures:
                    # optional argument provided so add as prefix
                    figname = "%s_%s.pdf" % (args.save_figures, key)
                else:
                    figname = "%s.pdf" % key
                self.figures[key].savefig(figname)

        if not args.no_show:
            plt.show()

    def _body(self, args):
        """This method contains the main body of code that defines the
        behaviour of the calculation. This is a simple example that
        invokes self.calculate and returns the result. The return
        value of this function is given to self.run which handles any
        outputs. Any figures created in this method should be added to
        the dictionary self.figures. The key used for each figure will
        be used as the basis of the name when saved. Any output tables
        created in this method should be added to the list self.tables.

        Parameters
        ----------
        args: argparse.Namespace object
            Namespace from argumentparser
        """
        return self.calculate(subset=(args.lower_bound, args.upper_bound))


def get_arg_parser():
    """Returns a generic argparser for all free energy calculation scripts"""
    parser = FEArgumentParser(add_help=False)
    parser.add_argument(
        '-d', '--directories', nargs='+', required=True, action='append',
        help="Location of folders containing ProtoMS output subdirectories. "
             "Multiple directories can be supplied to this flag and indicate "
             "repeats of the same calculation. This flag may be given "
             "multiple times and each instances is treated as an individual "
             "leg making up a single free energy difference e.g. vdw and ele "
             "contributions of a single topology calculation.")
    parser.add_argument(
        '--subdir', default='',
        help="Optional sub-directory for each lambda value to search within "
             "for simulation output. This is useful in, for instance, "
             "processing only the results of a GCAP calculation at a "
             "particular B value.")
    parser.add_argument(
        '-l', '--lower-bound', default=0., type=float,
        help="Value between 0 and 1 that determines the proportion to omit "
             "from the beginning of the simulation data series.")
    parser.add_argument(
        '-u', '--upper-bound', default=1., type=float,
        help="Value between 0 and 1 that determines the proportion to omit "
             "from the end of the simulation data series.")
    parser.add_argument(
        '-t', '--temperature', default=298.15, type=float,
        help='Temperature at which the simulation was run. Default=298.15K')
    parser.add_argument(
        '--pickle', help='Name of file in which to store results as a pickle.')
    parser.add_argument(
        '--save-figures', nargs='?', const='',
        help="Save figures produced by script. Takes optional argument that "
             "adds a prefix to figure names.")
    parser.add_argument(
        '--no-show', action='store_true', default=False,
        help="Do not display any figures on screen. Does not interfere with "
             "--save-figures.")
    return parser
