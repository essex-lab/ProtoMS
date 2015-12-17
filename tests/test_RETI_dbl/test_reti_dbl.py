"""This test is to check the ProtoMS setup of simulations with multiple lambda windows and the MPI implementation of protoms3 program."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# -------------------------------------------------------
# ProtoMS setup and RETI dual topology simulations test
# -------------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/RETI_dbl/"

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

output_files_setup = ["ethane_box.pdb", "eth-meo.tem", "run_free.cmd"]
output_subdirs = ["lam-0.000", "lam-0.330", "lam-0.670", "lam-1.000"]
out_sim_files = ["results", "accept", "all.pdb", "restart.prev", "warning", "info", "restart", "results_inst"]


class TestRETIdbl(unittest.TestCase):
    """Test for RETI/MPI function - short dual topology simulation."""
    def setUp(self):
        super(TestRETIdbl, self).setUp()

    def tearDown(self):
        super(TestRETIdbl, self).tearDown()

    def test_RETI_dbl(self):
        """Test for RETI/MPI function - short dual topology simulation."""

        if call("python2.7 $PROTOMSHOME/protoms.py -s dualtopology -l ethane.pdb methanol.pdb --nequil 0 --nprod 10 --lambdas 0.00 0.33 0.67 1.00 --ranseed 100000 --dumpfreq 1 --cleanup", shell=True) == 0:
            # Checking whether the required output files have been setup for RETI dual topology protoms.py setup.

            for out_files in output_files_setup:
                self.assertTrue(os.path.exists(out_files),
                                "ProtoMS setup output file {0} is missing.".format(out_files))

            # Checking content of RETI dual topology setup output files with reference data files.
            for out_files in output_files_setup:

                if out_files == "run_free.cmd":
                    if call("bash content_cmd_comp.sh", shell=True) == 0:
                        print("Content matched for command file {0}.".format(out_files))
                    else:
                        raise ValueError("Content mismatch between output and reference for file {0}".format(out_files))
                else:
                    self.assertTrue(filecmp.cmp(out_files, os.path.join(ref_dir, out_files)),
                                    "Content mismatch between output and reference for file {0}".format(out_files))

        else:
            raise simulationobjects.SetupError("ProtoMS RETI dual topology setup and command files generation failed")

        # Test for RETI dual topology simulation.
        if call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_free.cmd", shell=True) == 0:
            comparetools = testtools.CompareTools(ref_dir, verbose=True)
            # Checking whether the simulation output files have been created successfully for RETI dual topology.
            for subdir in output_subdirs:
                for outfile in out_sim_files:
                    outfile_rel = os.path.join("out_free", subdir, outfile)
                    self.assertTrue(os.path.exists(outfile_rel),
                                    "Simulation file {0} is missing.".format(outfile_rel))

                    # Checking content of RETI free phase leg of a dual topology simulation output files with reference data.
                    self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("RETI dual topology simulation was not successful.")

# Entry point to nosetests or unittests.
if __name__ == "__main__":
    unittest.main()
    nose.runmodule()
