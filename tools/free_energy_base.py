"""Collection of classes to form the basis of a replacement for the current
free energy calculation framework. """

from abc import abstractmethod
import argparse
import glob
import os
import numpy as np
import pymbar
from scipy.integrate import trapz
import simulationobjects as sim


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


class Estimator(object):
    """Base class for free energy estimators"""
    def __init__(self, lambdas):
        self.data = []
        self.lambdas = lambdas

    @abstractmethod
    def add_data(self):
        pass

    @abstractmethod
    def calculate(self):
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

    def subset(self, low_bound=0.0, high_bound=1.0, step=1):
        low_ind = int(len(self)*low_bound)
        high_ind = int(len(self)*high_bound)
        return self[low_ind:high_ind:step]


class TI(Estimator):
    """Estimate free energy differences using Thermodynamic Integration."""
    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.
        """
        self.data.append((series.forwfe-series.backfe) /
                         (series.lamf[0]-series.lamb[0]))
        return len(self.data[-1])

    def calculate(self, temp=300):
        """Calculate the free energy difference and return a PMF object."""
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
    """Top level class for performing a free energy calculation with
    simulation data.
    """
    def __init__(self, root_paths, estimators=[TI, BAR, MBAR]):
        """root_paths - a list of strings to ProtoMS output directories.
        The free energy will be calculated individually for each entry.
        estimators - a list of estimator classes to use"""
        self.root_paths = root_paths
        self.paths = [sorted(glob.glob(os.path.join(root_path, "lam-*")))
                      for root_path in self.root_paths]
        self.lambdas = [[float(path.split('lam-')[1]) for path in rep]
                        for rep in self.paths]
        self.estimators = {
            estimator: [estimator(lams)
                        for lams, rep in zip(self.lambdas, self.paths)]
            for estimator in estimators}
        self._extract_series(self.paths)

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
        """For each estimator return the evaluated potential of mean force."""
        return {estimator: [rep.subset(*subset).calculate() for rep in repeats]
                for estimator, repeats in self.estimators.iteritems()}

    def autoeqb(self):
        pass

    def test_equilibration(self, discard_limit):
        proportions = np.linspace(0.0, discard_limit, 10)
        return {prop: self.calculate(subset=(prop, 1.0))
                for prop in proportions}

    def test_convergence(self, discard_limit, lower_limit=0.):
        proportions = np.linspace(discard_limit, 1.0, 10)
        return {prop: self.calculate(subset=(lower_limit, prop))
                for prop in proportions}


class FEArgumentParser(argparse.ArgumentParser):
    """Thin wrapper around argparse.ArgumentParser designed to allow
    specification of individual arguments that cannot be used at the
    same time as one another."""
    def __init__(self, *args, **kwargs):
        self.clashes = {}
        try:
            for parent in kwargs['parents']:
                self.clashes.update(parent.clashes)
        except KeyError:
            pass
        argparse.ArgumentParser.__init__(self, *args, **kwargs)

    def parse_args(self, **kwargs):
        parsed = argparse.ArgumentParser.parse_args(self, **kwargs)
        self.check_clashes(parsed)
        return parsed

    def add_argument(self, *args, **kwargs):
        # figure out what the argument will be called in the parser namespace
        try:
            name = kwargs['dest']
        except KeyError:
            name = args[-1].strip('-')

        # store clashes if these are provided, otherwise empty list
        try:
            self.clashes[name] = kwargs.pop('clashes')
        except KeyError:
            self.clashes[name] = []

        argparse.ArgumentParser.add_argument(self, *args, **kwargs)

    def check_clashes(self, parsed):
        """Check that arguments in Namespace parsed are compatible,
        according to provided argument clashes."""
        for arg in self.clashes:
            # if parsed.arg has its default value it WAS NOT used so ignore
            if getattr(parsed, arg, None) == self.get_default(arg):
                continue
            for clash in self.clashes[arg]:
                # check clashes for this argument
                # if clashing argument has non-default value it WAS used
                # so object and quit
                if getattr(parsed, clash) != self.get_default(clash):
                    self.error('Cannot provide both --%s and --%s arguments' %
                               (arg, clash))


def get_arg_parser():
    """Returns a generic argparser for all free energy calculation scripts"""
    parser = FEArgumentParser(add_help=False)
    parser.add_argument('-d', '--directories', nargs='+', required=True,
                        help='output directories')
    # parser.add_argument(
    #     '--autoeqb', dest='autoeqb', action='store_true',
    #     help="Use automatic equilibration detection to determine how much "
    #          "data is included in the calculation")
    parser.add_argument(
        '-l', '--lower_bound', default=0., type=float,
        help="Define the lower bound of data to be used.")
    parser.add_argument(
        '-u', '--upper_bound', default=1., type=float,
        help="Define the upper bound of data to be used.")
    parser.add_argument(
        '--test_equilibration', default=None, type=float,
        help="Perform free energy calculations 10 times using varying "
             "proportions of the total data set provided. Data used will "
             "range from 100% of the dataset up to the proportion "
             "provided to this argument",
        clashes=('test_convergence', 'lower_bound', 'autoeqb'))
    parser.add_argument(
        '--test_convergence', default=None, type=float,
        help="Perform free energy calculations 10 times using varying "
             "proportions of the total data set provided. Data used will "
             "range from 100% of the dataset up to the proportion "
             "provided to this argument")
    return parser
