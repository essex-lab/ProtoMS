"""This test is to check the ProtoMS setup of simulations with multiple lambda windows and the MPI implementation of protoms3 program."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# --------------------------------------------------------
# ProtoMS setup and RETI single topology simulations test
# --------------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/RETI_sngl/"

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

output_files_setup = ["ethane_box.pdb", "ethtmeo_ele.tem", "ethtmeo_vdw.tem", "ethtmeo_comb.tem"]
outfiles_setup = ["run_comb_free.cmd", "run_comb_gas.cmd", "run_ele_free.cmd", "run_ele_gas.cmd", "run_vdw_free.cmd", "run_vdw_gas.cmd"]
output_subdirs = ["lam-0.000", "lam-0.330", "lam-0.670", "lam-1.000"]
out_sim_files = ["results", "accept", "all.pdb", "restart.prev", "warning", "restart", "results_inst"]


class TestRETIsngl(unittest.TestCase):
    """Test for RETI/MPI function."""

    def setUp(self):
        super(TestRETIsngl, self).setUp()

    def tearDown(self):
        super(TestRETIsngl, self).tearDown()

    def test_RETI_sngl(self):
        """Test for RETI single topology/MPI function."""
        comparetools = testtools.CompareTools(ref_dir, verbose=True)

        if call("python2.7 $PROTOMSHOME/protoms.py -s singletopology -l ethane.pdb methanol.pdb --nequil 0 --nprod 10 --lambdas 0.00 0.33 0.67 1.00 --ranseed 100000 --dumpfreq 1 --cleanup --singlemap single_cmap.dat", shell=True) == 0:
            # Checking whether the required output files have been setup for RETI single topology protoms.py setup.

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            for outfile in outfiles_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            # Checking content of RETI single topology setup output files with reference data files.
            for outfile in output_files_setup:
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

            for out_files in outfiles_setup:
                if call("bash content_cmd_comp.sh", shell=True) == 0:
                    print("Content matched for {0}.".format(out_files))
                else:
                    raise ValueError("Content mismatch between output and reference for file {0}".format(out_files))

        else:
            raise simulationobjects.SetupError("ProtoMS setup for RETI single topology and command files generation failed.")

        if call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_comb_free.cmd", shell=True) == 0:
            # Checking whether the simulation output files have been created successfully for RETI free phase leg of a single topology for combined perturbation.
            for subdir in output_subdirs:
                for outfile in out_sim_files:
                    outfile_rel = os.path.join("out_comb_free", subdir, outfile)
                    self.assertTrue(os.path.exists(outfile_rel),
                                    "Simulation file {0} is missing.".format(outfile_rel))

                    # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                    self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("RETI free phase leg of a single topology simulation is not successful.")

        if call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_comb_gas.cmd", shell=True) == 0:
            # Checking whether the simulation output files have been created successfully for RETI gas phase leg of a single topology for combined perturbation.
            for subdir in output_subdirs:
                for outfile in out_sim_files:
                    outfile_rel = os.path.join("out_comb_free", subdir, outfile)
                    self.assertTrue(os.path.exists(outfile_rel),
                                    "Simulation file {0} is missing.".format(outfile_rel))

                    # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                    self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("RETI gas phase leg of a single topology simulation was not successful.")


# Entry point to nosetests or unittests
if __name__ == "__main__":
    unittest.main()
    nose.runmodule()
