import pymbar
import numpy as np
import unittest
from free_energy_base import BAR, MBAR, TI
from calc_ti_decomposed import TI_decomposed
from simulationobjects import boltz, SnapshotResults


class TestBAR(unittest.TestCase):
    estimator_class = BAR

    def setUp(self):
        """Create a test based on the pymbar HarmonicOscillatorsTestCase.
        """
        self.test = pymbar.testsystems.HarmonicOscillatorsTestCase(
            O_k=[0, 0, 0, 0, 0],
            K_k=[4., 5., 6., 7., 8.])
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

        self.estimator = self.estimator_class(list(self.lambdas))
        self.add_data()

    def add_data(self):
        for series in self.series:
            self.estimator.add_data(series)

    def test_free_energies(self):
        """Check that free energies returned by the estimator's calculate
        method match the analytical results for the system.
        """

        # default test parameters use thermodynamic beta = 1
        pmf = self.estimator.calculate(temp=1/boltz)
        for fe1, fe2 in zip(pmf.values, self.test.analytical_free_energies()):
            self.assertAlmostEqual(fe1, fe2, places=1)

    def test_getitem(self):
        lens = [len(series)
                for series in self.estimator._get_data()]

        for dN in (10, 20, 50, 100, 500):
            new_est = self.estimator[dN:]
            sliced_lens = [len(series)
                           for series in new_est._get_data()]
            print lens, sliced_lens
            for l1, l2 in zip(lens, sliced_lens):
                self.assertEqual(l1, l2 + dN)

        # check that new_est is a different instance
        self.assertIsNot(self.estimator, new_est)


class TestMBAR(TestBAR):
    estimator_class = MBAR

    def test_getitem(self):
        lens = [series.shape[1]
                for series in self.estimator._get_data()]

        for dN in (10, 20, 50, 100, 500):
            new_est = self.estimator[dN:]
            sliced_lens = [series.shape[1]
                           for series in new_est._get_data()]
            print lens, sliced_lens
            for l1, l2 in zip(lens, sliced_lens):
                self.assertEqual(l1, l2 + dN)

        # check that new_est is a different instance
        self.assertIsNot(self.estimator, new_est)


class TestTI(TestBAR):
    estimator_class = TI

    def add_data(self):
        """For the Harmonic oscillator it is simple to work out the gradient
        of the energy with respect to lambda under certain
        conditions. This data is spoofed into series.forwfe such that
        when the estimator calculates the gradient it recovers the
        correct values. Requires strictly linear scaling of spring
        constants, and for all potentials to have the same bond length!

        dU(\lambda)/d\lambda = (k_1-k_0)/2*(x-b_0)^2
        """
        dlam = 0.001
        for series, state in zip(self.series, self.u_kln):
            series.lamb = [series.lam[0]-dlam]
            series.lamf = [series.lam[0]+dlam]
            series.forwfe = dlam*2*(series.feenergies[1.0] -
                                    series.feenergies[0.0])
            series.backfe = np.zeros_like(series.forwfe)
            self.estimator.add_data(series)
