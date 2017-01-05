"""This test is to check the protoms.py setup and protoms3 run of gcmc simulations."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# ---------------------------------------------
# ProtoMS setup and GCMC simulations test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = os.path.join(proto_env, "tests/gcmc/")

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

out_gcmc_tools = ["gcmc_box.pdb"]
output_files_setup = ["gcmc_wat.pdb", "run_gcmc.cmd"]
output_subdirs = ["b_+19.000","b_+20.000"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning"]


class TestGCMC(unittest.TestCase):
    """Test for GCMC function."""

    def setUp(self):
        super(TestGCMC, self).setUp()

    def tearDown(self):
        super(TestGCMC, self).tearDown()

    def test_gcmc(self):
        """Test for GCMC function."""

        if((call("python2.7 $PROTOMSHOME/tools/make_gcmcbox.py -s wat.pdb", shell=True)) == 0):
            # Checking whether GCMC tool created the box successfully.

            for outfile in out_gcmc_tools:
                self.assertTrue(os.path.exists(outfile),
                                "GCMC setup output file {0} is missing.".format(outfile))

            for outfile in out_gcmc_tools:
                # Checking content of GCMC Box with reference GCMC Box file.
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

        else:
            raise simulationobjects.SetupError("Creation of GCMC simulation box was not successful!")

        if((call("python2.7 $PROTOMSHOME/protoms.py -sc protein.pdb -s gcmc --gcmcwater wat.pdb --gcmcbox gcmc_box.pdb --adams 19 20 --nequil 0 --nprod 1000 --ranseed 100000 --dumpfreq 100 --capradius 26 -w water.pdb", shell=True)) == 0):
            # Checking whether the required output files have been setup for GCMC simulation.

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "Setup output file {0} is missing.".format(outfile))

            for outfile in output_files_setup:
                # Checking content of GCMC setup output files with reference data in files.
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

            print("Setup and command files generation is successful.")
        else:
            raise simulationobjects.SetupError("ProtoMS setup and command files generation is not successful.")

        # Test for ProtoMS simulation.
        if((call("mpirun -np 2 $PROTOMSHOME/build/protoms3 run_gcmc.cmd", shell=True)) == 0):
            # Checking whether the simulation output files have been created successfully.
            comparetools = testtools.CompareTools(ref_dir, verbose=True)
            for subdir in output_subdirs:
                for outfile in out_sim_files:
                    outfile_rel = os.path.join("out_gcmc", subdir, outfile)
                    self.assertTrue(os.path.exists(outfile_rel),
                                    "GCMC simulation file {0} is missing.".format(outfile_rel))

                    # Checking content of GCMC simulation output files with reference data in files.
                    self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("GCMC simulation was not successful.")


# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
