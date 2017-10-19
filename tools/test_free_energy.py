import pymbar
import numpy as np
import unittest
from free_energy_base import BAR
from simulationobjects import boltz, SnapshotResults

beta = 1. / (boltz*300)


class TestEstimators(unittest.TestCase):

    def setUp(self):
        self.test = pymbar.testsystems.HarmonicOscillatorsTestCase()
        # k = sample state index
        # l = eval state index
        # n = sample index
        self.x_kn, self.u_kln, self.N_k = self.test.sample(N_k=[10E4]*5,
                                                           mode='u_kln')

        self.lambdas = np.linspace(0., 1., len(self.N_k))
        self.series = []
        for lam, state in zip(self.lambdas, self.u_kln):
            snapshot = SnapshotResults()
            snapshot.feenergies = {lam: state[i]
                                   for i, lam in enumerate(self.lambdas)}
            snapshot.lam = [lam]
            self.series.append(snapshot)

    def testBAR(self):

        est = BAR(list(self.lambdas))

        for series in self.series:
            est.add_data(series)

        pmf = est.calculate(temp=1/boltz)
        for fe1, fe2 in zip(pmf.values, self.test.analytical_free_energies()):
            self.assertAlmostEqual(fe1, fe2, places=1)


