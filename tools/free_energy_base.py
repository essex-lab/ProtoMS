"""Collection of classes to form the basis of a replacement for the current
free energy calculation framework. """

from abc import abstractmethod, ABCMeta
from copy import copy
import glob
from operator import add
import os
import matplotlib
import numpy as np
import pickle
import pymbar
from scipy.integrate import trapz
import simulationobjects as sim
from free_energy_argument_parser import FEArgumentParser

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['lines.linewidth'] = 3


class PMF(object):
    """A potential of mean force, describing the free energy change as a
    function of the lambda value.
    """
    def __init__(self, lambdas, values):
        """lambdas - a sequence containing the lambda values at which the free
        energy has been evaluated.
        values - sequence of the free energy values at each lambda
        """
        self.lambdas = lambdas
        self.values = values

    @property
    def dG(self):
        """Return the free energy difference at the PMF end points."""
        return self.values[-1] - self.values[0]

    def __neg__(self):
        return PMF(self.lambdas, [-val for val in self.values])

    def __iter__(self):
        return iter(self.values)


class CompositePMF(PMF):
    def __init__(self, *args):
        if len({tuple(arg.lambdas) for arg in args}) > 1:
            raise ValueError("All data must use the same lambda values")
        self.lambdas = args[0].lambdas
        self.values = []
        for dat in zip(*args):
            try:
                # if constructing from PMFs just have floats
                self.values.append(FreeEnergy(*_get_mean_std_err(dat)))
            except TypeError:
                # constructing from CompositePMFs already have FreeEnergy's
                self.values.append(np.sum(dat))

    def __add__(self, other):
        return CompositePMF(self, other)

    def __neg__(self):
        other = copy(self)
        other.values = [-val for val in other.values]
        return other

    def plot(self, ax, **kwargs):
        y = np.array([fe.value for fe in self.values])
        err = np.array([fe.error for fe in self.values])

        line = ax.plot(self.lambdas, y, **kwargs)[0]
        ax.plot(self.lambdas, y+err, '--', linewidth=1, color=line.get_color())
        ax.plot(self.lambdas, y-err, '--', linewidth=1, color=line.get_color())


class Result(object):
    def __init__(self, *args):
        self.data = args
        if len({tuple(pmf.lambdas) for dat in self.data for pmf in dat}) > 1:
            raise ValueError("All data must use the same lambda values")

    @property
    def lambdas(self):
        try:
            return self.data[0][0].lambdas
        except IndexError:
            raise ValueError("This result does not contain any data.")

    @property
    def dG(self):
        """Return a free energy object from this result's data."""
        data_dGs = [[pmf.dG for pmf in dat] for dat in self.data]
        FEs = [FreeEnergy(*_get_mean_std_err(dG)) for dG in data_dGs]
        return sum(FEs, FreeEnergy(0., 0.))

    @property
    def pmf(self):
        pmfs = [CompositePMF(*dat) for dat in self.data]
        return reduce(add, pmfs)

    def __add__(self, other):
        return Result(*(self.data + other.data))

    def __sub__(self, other):
        return self + -other

    def __neg__(self):
        return Result(*[[-pmf for pmf in dat] for dat in self.data])


class FreeEnergy(object):
    def __init__(self, value, error):
        self.value = value
        self.error = error

    def __add__(self, other):
        return FreeEnergy(self.value + other.value,
                          (self.error**2 + other.error**2)**0.5)

    def __sub__(self, other):
        return FreeEnergy(self.value - other.value,
                          (self.error**2 + other.error**2)**0.5)

    def __neg__(self):
        return FreeEnergy(-self.value, self.error)

    def __eq__(self, other):
        return (self.value == other.value) and (self.error == other.error)

    def __str__(self):
        return "%9.4f +/- %.4f" % (self.value, self.error)

    def __repr__(self):
        return "<FreeEnergy: value=%.4f error=%.4f>" % (self.value, self.error)


class Estimator(object):
    __metaclass__ = ABCMeta
    """Base class for free energy estimators"""
    def __init__(self, lambdas):
        self.data = []
        self.lambdas = lambdas

    @abstractmethod
    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.
        """
        pass

    @abstractmethod
    def calculate(self, temp=300):
        """Calculate the free energy difference and return a PMF object.

        temp: float, optional
              temperature of calculation
        """
        pass

    def __getitem__(self, val):
        """Return a new class instance with series[val] applied to each
        individual data series.
        """
        new_est = self.__class__(self.lambdas)
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
        """Return a new class instance containing truncated data series."""
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
        self.data.append((series.forwfe-series.backfe) /
                         (series.lamf[0]-series.lamb[0]))
        return len(self.data[-1])

    def calculate(self, temp=300):
        gradients = [gradient_data.mean() for gradient_data in self.data]
        pmf_values = [trapz(gradients[:i], self.lambdas[:i])
                      for i in xrange(1, len(self.lambdas) + 1)]
        return PMF(self.lambdas, pmf_values)


class BAR(Estimator):
    """Estimate free energy differences using Bennett's Acceptance Ratio."""
    def add_data(self, series):
        lam = series.lam[0]
        lam_ind = self.lambdas.index(lam)
        lamf = self.lambdas[lam_ind+1] if lam != 1.0 else 1.0
        lamb = self.lambdas[lam_ind-1] if lam != 0.0 else 0.0
        self.data.append(
            np.array([series.feenergies[lam] - series.feenergies[lamb],
                      series.feenergies[lam] - series.feenergies[lamf]]))
        return len(self.data[-1][0])

    def calculate(self, temp=300):
        """Calculate the free energy difference and return a PMF object."""
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
        self.data.append(np.array([series.feenergies[lam]
                                   for lam in sorted(series.feenergies)]))
        return len(self.data[-1][0])

    def calculate(self, temp=300):
        beta = 1./(sim.boltz*temp)
        mbar = pymbar.MBAR(np.array(self.data)*beta,
                           [len(dat[0]) for dat in self.data])
        FEs = mbar.getFreeEnergyDifferences()[0]/beta
        return PMF(self.lambdas,
                   [FEs[0, i] for i in xrange(len(self.data))])


class FreeEnergyCalculation(object):
    """Top level class for performing a free energy calculation from
    ProtoMS simulation outputs.
    """

    def __init__(self, root_paths, temperature, estimators=[TI, BAR, MBAR]):
        """root_paths - a list of strings to ProtoMS output directories.
        The free energy will be calculated individually for each entry.
        temperature - simulation temperature in Kelvin.
        estimators - a list of estimator classes to use."""
        self.root_paths = root_paths
        self.temperature = temperature
        self.paths = [sorted(glob.glob(os.path.join(root_path, "lam-*")))
                      for root_path in self.root_paths]
        self.lambdas = [[float(path.split('lam-')[1]) for path in rep]
                        for rep in self.paths]
        self.estimators = {
            estimator: [estimator(lams)
                        for lams, rep in zip(self.lambdas, self.paths)]
            for estimator in estimators}
        self._extract_series(self.paths)
        self.figures = {}

    def _extract_series(self, paths):
        """For all results files extract the data series and supply these
        to each estimator instance."""
        for i, repeat in enumerate(paths):
            min_len = 10E20
            for path in repeat:
                rf = sim.ResultsFile()
                rf.read(os.path.join(path, 'results'))
                rf.make_series()
                for cls in self.estimators:
                    data_len = self.estimators[cls][i].add_data(rf.series)
                    min_len = data_len if data_len < min_len else min_len

            # in cases where calculations terminate prematurely there
            # can be slight differences in the length of data series
            # within a repeat. Here we standardise the length for
            # later convenience
            for cls in self.estimators:
                self.estimators[cls][i] = self.estimators[cls][i][:min_len]

    def calculate(self, subset=(0., 1., 1)):
        """For each estimator return the evaluated potential of mean force.

        Returns
        -------
        dict:
          results of calculation stored in the below structure:
          return_val[estimator_class][i] = PMF instance
          where i is an list index indicating a repeat
        """
        results = {}
        for est, repeats in self.estimators.items():
            results[est] = Result(
                [rep.subset(*subset).calculate(self.temperature)
                 for rep in repeats])
        return results

    def run(self, args):
        results = self._body(args)

        if args.pickle is not None:
            with open(args.pickle, 'w') as f:
                pickle.dump(results, f, protocol=2)

        # slightly tortured logic arises below from the behaviour of argparse
        if (args.save_figures is None) or args.save_figures:
            for key in self.figures:
                if args.save_figures is None:
                    figname = "%s.pdf" % key
                else:
                    figname = "%s_%s.pdf" % (args.save_figures, key)
                self.figures[key].savefig(figname)

        plt.show()

    def _body(self, args):
        """This method contains the main body of code that defines the
        behaviour of the calculation. This is a simple example that
        invokes self.calculate and returns the result. The return
        value of this function is given to self.run which handles any
        outputs. Any figures created in this method should be added to
        the dictionary self.figures. The key used for each figure will
        be used as the basis of the name when saved.
        """
        return self.calculate(subset=(args.lower_bound, args.upper_bound))


def _get_mean_std_err(dat):
    return np.mean(dat), np.std(dat)/len(dat)**0.5


def get_arg_parser():
    """Returns a generic argparser for all free energy calculation scripts"""
    parser = FEArgumentParser(add_help=False)
    parser.add_argument('-d', '--directories', nargs='+', required=True,
                        help='output directories')
    parser.add_argument(
        '-l', '--lower-bound', default=0., type=float,
        help="Define the lower bound of data to be used.")
    parser.add_argument(
        '-u', '--upper-bound', default=1., type=float,
        help="Define the upper bound of data to be used.")
    parser.add_argument(
        '-t', '--temperature', default=298.15, type=float,
        help='Temperature at which the simulation was run. Default=298.15K')
    parser.add_argument(
        '--pickle', help='Name of file in which to store results as a pickle.')
    parser.add_argument(
        '--save-figures', nargs='?', default='',
        help="Save figures produced by script. Takes optional argument that "
             "adds a prefix to figure names.")
    return parser
