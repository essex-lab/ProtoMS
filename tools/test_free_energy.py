import pymbar
import numpy as np
import unittest
from free_energy_base import BAR, MBAR, TI
from simulationobjects import boltz, SnapshotResults


class TestBAR(unittest.TestCase):
    estimator_class = BAR

    def setUp(self):
        """Create a test based on the pymbar HarmonicOscillatorsTestCase.
        """
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

    def test_slicing(self):
        lens = [len(series)
                for dat in self.estimator.data for series in dat]
        dN = 50
        slc = slice(dN, None, None)
        self.estimator.apply_slice(slc)
        sliced_lens = [len(series)
                       for dat in self.estimator.data for series in dat]
        
        for l1, l2 in zip(lens, sliced_lens):
            self.assertEqual(l1, l2 + dN)
        
    # def test_free_energies(self):
    #     est = BAR(list(self.lambdas))
    #     self.BARs(est)


class TestMBAR(TestBAR):
    estimator_class = MBAR
    # def test_free_energies(self):
    #     est = MBAR(list(self.lambdas))
    #     self.BARs(est)

    # def test_slicing(self):
    #     lens = [len(series)
    #             for dat in self.estimator.data for series in dat]
    #     dN = 50
    #     slc = slice(dN, None, None)
    #     self.estimator.apply_slice(slc)
    #     sliced_lens = [len(series)
    #                    for dat in self.estimator.data for series in dat]
        
    #     for l1, l2 in zip(lens, sliced_lens):
    #         self.assertEqual(l1, l2 + dN)


class TestTI(TestBAR):
    estimator_class = TI

    def add_data(self):
        """For the time being the data added is meaningless but has the same
        form and hence is suitable for testing slicing"""
        for series in self.series:
            self.estimator.data.append(series.feenergies[0.0])
        # for series in self.series:
        #     self.estimator.add_data(series)

    def test_slicing(self):
        lens = [len(series)
                for series in self.estimator.data]
        dN = 50
        slc = slice(dN, None, None)
        self.estimator.apply_slice(slc)
        sliced_lens = [len(series)
                       for series in self.estimator.data]

        for l1, l2 in zip(lens, sliced_lens):
            self.assertEqual(l1, l2 + dN)

    def test_free_energies(self):
        pass

