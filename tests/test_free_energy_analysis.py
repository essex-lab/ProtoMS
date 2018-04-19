"""Execute with 'python2.7 test_free_energy.py'"""

import pymbar
import numpy as np
import os
import shutil
import sys
import unittest
import framework

sys.path.append(os.path.join(os.environ['PROTOMSHOME'], 'tools'))

import calc_dg2
import calc_ti_decomposed
import calc_dg_cycle
import calc_gci2
from free_energy_base import BAR, MBAR, TI, get_alchemical_arg_parser
from free_energy_base import PMF, Quantity, Result
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
            self.assertAlmostEqual(fe1.value, fe2, places=1)

    def test_getitem(self):
        lens = [series.shape[-1]
                for series in self.estimator.data]

        for dN in (10, 20, 50, 100, 500):
            new_est = self.estimator[dN:]
            sliced_lens = [series.shape[-1]
                           for series in new_est.data]

            for l1, l2 in zip(lens, sliced_lens):
                self.assertEqual(l1, l2 + dN)

        # check that new_est is a different instance
        self.assertIsNot(self.estimator, new_est)


class TestMBAR(TestBAR):
    estimator_class = MBAR


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


class TestArgumentParsers(unittest.TestCase):
    """Check for errors in the inheritance scheme used for argument parsers"""
    def test_base_parser(self):
        get_alchemical_arg_parser()

    def test_dg2_parser(self):
        calc_dg2.get_arg_parser()

    def test_ti_decomposed(self):
        calc_ti_decomposed.get_arg_parser()

    def test_dg_cycle(self):
        calc_dg_cycle.get_arg_parser()


class TestPMF(unittest.TestCase):
    def setUp(self):
        self.pmf_values = np.linspace(0., 5., 16)
        self.lambdas = np.linspace(0., 1., 16)
        self.pmfs = [list(self.pmf_values*i*0.05) for i in xrange(1, 3)]
        self.single = PMF(self.lambdas, self.pmfs[0])
        self.double = PMF(self.lambdas, *self.pmfs)

    def test_dG(self):
        self.assertEqual(self.single.dG,
                         Quantity(self.pmf_values[-1]*0.05, 0.))
        self.assertEqual(self.double.dG.value, self.pmf_values[-1]*0.05*1.5)

    def test_negation(self):
        """Test use of negation operator"""
        self.single = PMF(self.lambdas, self.pmfs[0])
        self.double = PMF(self.lambdas, *self.pmfs)

        # as copy operator has been used check that data of original
        # object is not changed by the copied version
        pre = self.single.dG
        -self.single
        self.assertNotEqual(pre, -self.single.dG)

        # check for expected behaviours
        self.assertEqual((-self.single).dG, -(self.single.dG))
        self.assertEqual((-self.double).dG, -(self.double.dG))

    def test_addition(self):
        """Test addition"""
        self.single = PMF(self.lambdas, self.pmfs[0])
        self.double = PMF(self.lambdas, *self.pmfs)

        self.assertEqual((self.single+self.single).dG.value,
                         self.single.dG.value*2)
        self.assertEqual((self.double+self.single).dG.value,
                         self.single.dG.value*2.5)

        # test lambda value checking on addition
        self.single.coordinate = [0, 1]
        with self.assertRaises(ValueError):
            self.single + self.double


class TestResult(unittest.TestCase):
    def setUp(self):
        self.pmf_values = np.linspace(0., 5., 16)
        lambdas = np.linspace(0., 1., 16)
        self.pmfs = [PMF(lambdas, list(self.pmf_values*i*0.05))
                     for i in xrange(1, 4)]
        self.pmfs2 = [PMF(lambdas, list(self.pmf_values*2*i*0.05))
                      for i in xrange(1, 4)]
        self.single_result = Result(self.pmfs)
        self.double_result = Result(self.pmfs, self.pmfs2)

    def test_init(self):
        """Test that validation of data lambda values during construction
        works correctly. Must all be the same else a ValueError should be
        raised."""
        self.pmfs[0].coordinate = np.linspace(0., 2., 16)
        with self.assertRaises(ValueError):
            Result(self.pmfs)
        with self.assertRaises(ValueError):
            Result(self.pmfs, self.pmfs2)

    def test_dG(self):
        """Test that expected free energy values are returned."""
        self.assertEqual(self.single_result.dG.value,
                         self.pmf_values[-1]*2*0.05)
        self.assertEqual(self.double_result.dG.value,
                         self.pmf_values[-1]*6*0.05)

    def test_pmf(self):
        self.assertEqual(
            [round(fe.value, 8) for fe in self.single_result.pmf],
            [round(val*2*0.05, 8) for val in self.pmf_values])
        self.assertEqual(
            [round(fe.value, 8) for fe in self.double_result.pmf],
            [round(val*6*0.05, 8) for val in self.pmf_values])

    def test_negation(self):
        """Test use of negation operator"""
        self.assertEqual((-self.single_result).dG, -(self.single_result.dG))
        self.assertEqual((-self.double_result).dG, -(self.double_result.dG))

    def test_addition(self):
        """Test addition of Results objects"""
        self.assertEqual((self.single_result+self.single_result).dG.value,
                         self.single_result.dG.value*2)
        self.assertEqual((self.double_result+self.single_result).dG.value,
                         self.single_result.dG.value*4)


class testCalcDg(framework.BaseTest):
    """Test main driver classes associated with each of the analysis
    scripts.

    As very artificial data is used to construct these tests we are mainly
    looking for them to pass without throwing exceptions."""
    ref_dir = "tests/free_energy_analysis/"

    input_files = ['ala_gly']
    executable = ''
    args = []
    output_files = ['results.pkl', 'pref_pmf.pdf']
    cmdline = (
        "-d ala_gly/out1_free ala_gly/out2_free -d ala_gly/out3_free --pmf"
        " --no-show --pickle results.pkl --save-figures pref"
    )

    @classmethod
    def _helper_copy_input_files(self):
        for filename in self.input_files:
            try:
                shutil.copytree(os.path.join(self._full_ref_dir, filename),
                                filename)
            except IOError:
                raise IOError(
                    "The required reference input file {0} could not be "
                    "copied".format(filename))

    @classmethod
    def _helper_clean_files(cls):
        for filename in cls.input_files + cls.output_files:
            try:
                os.remove(filename)
            except OSError:
                pass
            try:
                shutil.rmtree(filename)
            except OSError:
                pass

    def tearDown(self):
        for filename in self.output_files:
            self.assertTrue(
                os.path.exists(filename),
                "Expected output file {0} is missing".format(filename))

    def test(self):
        calc_dg2.run_script(self.cmdline.split())


class testCalcDgCycleDual(testCalcDg):
    input_files = ['ala_gly', 'gly_val', 'ala_val']
    output_files = ['results.pkl']
    cmdline = ("-d ala_gly gly_val ala_val --signs + + - --dual"
               " --no-show --pickle results.pkl")

    def test(self):
        calc_dg_cycle.run_script(self.cmdline.split())


class testCalcDgCycleComb(testCalcDgCycleDual):
    cmdline = ("-d ala_gly gly_val ala_val --signs + + - "
               "--single comb --no-show --pickle results.pkl")


class testCalcDgCycleSep(testCalcDgCycleDual):
    cmdline = ("-d ala_gly gly_val ala_val --signs + + - "
               "--single sep --no-show --pickle results.pkl")


class testCalcTiDecomposed(testCalcDg):
    input_files = ['ala_gly']
    output_files = ['results.pkl', 'pref_decomposed_pmfs.pdf',
                    'pref_decomposed.pdf']
    cmdline = ("-d ala_gly/out1_free ala_gly/out2_free --pmf --no-show "
               "--pickle results.pkl --save-figures pref")

    def test(self):
        calc_ti_decomposed.run_script(self.cmdline.split())


class testCalcTiDecomposedDual(testCalcTiDecomposed):
    cmdline = ("-d ala_gly/out1_free ala_gly/out2_free --dual --pmf "
               " --no-show --pickle results.pkl --save-figures pref")


class testCalcDgGCAP(testCalcDg):
    input_files = ['gcap']
    output_files = ['results.pkl', 'pref_pmf_2d.pdf', 'pref_pmf.pdf']
    cmdline = (
        "-d gcap/ --est gcap bar --pmf --subdir b_-9.700 --name"
        " results_inst --no-show --pickle results.pkl --save-figures pref"
    )


class testCalcGCI(testCalcDg):
    input_files = ['gcap/lam-1.000']
    output_files = ['results.pkl']
    cmdline = ("-d gcap/lam-1.000 --save-figures pref --no-show"
               " --pickle results.pkl --name results_inst -v 30.")

    def test(self):
        calc_gci2.run_script(self.cmdline.split())


class testCalcGCIZero(testCalcGCI):
    """This test works with data that has a uniform water occupancy of zero.
    This ensures the corner case of a flat titration curve is
    handled gracefully."""
    input_files = ['gcap/lam-0.000']
    output_files = ['results.pkl']
    cmdline = ("-d gcap/lam-0.000 --save-figures pref --no-show"
               " --pickle results.pkl --name results_inst -v 30.")


if __name__ == '__main__':
    unittest.main()
