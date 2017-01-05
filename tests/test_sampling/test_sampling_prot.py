"""This test is to check the Sampling MC moves."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# ---------------------------------------------
# ProtoMS Equilibration Sampling MC moves test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/sampling/"

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

output_files_setup = ["run_bnd.cmd"]
out_sim_files = ["all.pdb", "warning", "results", "restart.prev", "restart", "accept"]


class TestSampling(unittest.TestCase):

    """Test for Sampling function."""

    def setUp(self):
        super(TestSampling, self).setUp()

    def tearDown(self):
        super(TestSampling, self).tearDown()

    def test_sampling(self):
        """ Test for Sampling function."""
        comparetools = testtools.CompareTools(ref_dir, verbose=True)

        if call("python2.7 $PROTOMSHOME/protoms.py -s sampling -l dcb.pdb -p protein.pdb --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 --gaff gaff14", shell=True) == 0:

            # Checking whether the required output files have been setup for Sampling MC moves.
            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "Setup output file {0} is missing.".format(outfile))

            for outfile in output_files_setup:
                # Checking content of setup output files with reference data in files.
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

        else:
            raise simulationobjects.SetupError("ProtoMS setup for sampling MC moves is not successful!")

        if call("$PROTOMSHOME/build/protoms3 run_bnd.cmd", shell=True) == 0:
            # Checking whether the simulation output files have been created successfully for Sampling MC moves.
            for outfile in out_sim_files:
                outfile_rel = os.path.join("out_bnd", outfile)
                self.assertTrue(os.path.exists(outfile_rel),
                                "Simulation file {0} is missing.".format(outfile_rel))

                # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("Sampling MC protoms simulation failed!")

# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
