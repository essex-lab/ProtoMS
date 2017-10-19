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
    function of the lambda value."""
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

    @abstractmethod
    def apply_slice(self):
        pass


class TI(Estimator):
    """Estimate free energy differences using Thermodynamic Integration."""
    def add_data(self, series):
        """Save data from a SnapshotResults.series object
        for later calculation"""
        self.data.append((series.forwfe-series.backfe) /
                         (series.lamf[0]-series.lamb[0]))

    def calculate(self, temp=300):
        """Calculate the free energy difference and return a PMF object."""
        gradients = [gradient_data.mean() for gradient_data in self.data]
        pmf_values = [trapz(gradients[:i], self.lambdas[:i])
                      for i in xrange(1, len(self.lambdas) + 1)]
        return PMF(self.lambdas, pmf_values)

    def apply_slice(self, slc):
        """Apply provided slice object to the estimators data series."""
        for i, series in enumerate(self.data):
            self.data[i] = series[slc]


class BAR(Estimator):
    """Estimate free energy differences using Bennett's Acceptance Ratio."""
    def add_data(self, series):
        lam = series.lam[0]
        lam_ind = self.lambdas.index(lam)
        lamf = self.lambdas[lam_ind+1] if lam != 1.0 else 1.0
        lamb = self.lambdas[lam_ind-1] if lam != 0.0 else 0.0
        self.data.append([series.feenergies[lam] - series.feenergies[lamb],
                          series.feenergies[lam] - series.feenergies[lamf]])

    def calculate(self, temp=300):
        beta = 1./(sim.boltz*temp)
        pmf_values = [0.0]
        for low_lam, high_lam in zip(self.data, self.data[1:]):
            pmf_values.append(pmf_values[-1] +
                              pymbar.BAR(-low_lam[1]*beta,
                                         -high_lam[0]*beta)[0]/beta)
        return PMF(self.lambdas, pmf_values)

    def apply_slice(self, slc):
        for i, dat in enumerate(self.data):
            for j, series in enumerate(dat):
                self.data[i][j] = series[slc]


class MBAR(Estimator):
    """Estimate free energy differences using the Multistate Bennett's
    Acceptante Ratio.
    """
    def add_data(self, series):
        self.data.append(np.array([series.feenergies[lam]
                                   for lam in sorted(series.feenergies)]))

    def calculate(self, temp=300):
        beta = 1./(sim.boltz*temp)
        mbar = pymbar.MBAR(np.array(self.data)*beta,
                           [len(dat[0]) for dat in self.data])
        FEs = mbar.getFreeEnergyDifferences()[0]/beta
        return PMF(self.lambdas,
                   [FEs[0, i] for i in xrange(len(self.data))])

    def apply_slice(self, slc):
        for i, dat in enumerate(self.data):
            self.data[i] = dat[:, slc]


class FreeEnergyCalculation(object):
    """Top level class for performing a free energy calculation with
    simulation data."""
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
        self._extract_ser(self.paths)

    def _extract_ser(self, paths):
        """For all results files extract the data series and supply these
        to each estimator instance."""
        for i, repeat in enumerate(paths):
            for path in repeat:
                rf = sim.ResultsFile()
                rf.read(os.path.join(path, 'results'))
                rf.make_series()
                for estimator in self.estimators:
                    self.estimators[estimator][i].add_data(rf.series)

    def calculate(self):
        """For each estimator return the evaluated potential of mean force."""
        return {estimator: [rep.calculate() for rep in repeats]
                for estimator, repeats in self.estimators.iteritems()}


def get_arg_parser():
    """Returns a generic argparser for all free energy calculation scripts"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directories', nargs='+',
                        help='output directories')
    return parser
