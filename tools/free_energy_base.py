"""Collection of classes to form the basis of a replacement for the current
free energy calculation framework. """

import glob
import os
from scipy.integrate import trapz
import simulationobjects as sim
import numpy as np


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
    def __init__(self):
        self.data = []


class TI(Estimator):
    """Estimate free energy differences using Thermodynamic Integration."""
    def add_data(self, series):
        """Save data from a SnapshotResults.series object
        for later calculation"""
        self.data.append((series.forwfe-series.backfe) /
                         (series.lamf[0]-series.lamb[0]))

    def calculate(self, lambdas):
        """Calculate the free energy difference and return a PMF object."""
        gradients = [gradient_data.mean() for gradient_data in self.data]
        pmf_values = [trapz(gradients[:i], lambdas[:i])
                      for i in xrange(1, len(lambdas) + 1)]
        return PMF(lambdas, pmf_values)


class FreeEnergyCalculation(object):
    """Top level class for performing a free energy calculation with
    simulation data."""
    def __init__(self, root_paths, estimators=[TI]):
        """root_paths - a list of strings to ProtoMS output directories. 
        The free energy will be calculated individually for each entry.
        estimators - a list of estimator classes to use"""
        self.root_paths = root_paths
        self.paths = [sorted(glob.glob(os.path.join(root_path, "lam-*")))
                      for root_path in self.root_paths]
        self.lambdas = [[float(path.split('lam-')[1]) for path in rep]
                        for rep in self.paths]
        self.estimators = {estimator: [estimator() for rep in self.paths]
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
        return {estimator: [rep.calculate(self.lambdas[0])
                            for rep, lambdas in zip(repeats, self.lambdas)]
                for estimator, repeats in self.estimators.iteritems()}


# calc = FreeEnergyCalculation(['/tmp/out_comb_bnd'],
#                              estimators=[TI, TI_decomposed])
# results = calc.calculate()
