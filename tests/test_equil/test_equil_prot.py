"""This test is to check the equilibration MC moves prior to Sampling MC moves."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# ---------------------------------------------
# ProtoMS Equilibration MC moves test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

ref_dir = proto_env + "/tests/equil/"
input_files_setup = ["dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb"]
output_files_setup = ["run_bnd.cmd"]
out_sim_files = ["equil_bnd.pdb", "warning"]


class TestEquilibrationSetup(unittest.TestCase):
    """Test for equilibration function."""

    def setUp(self):
        super(TestEquilibrationSetup, self).setUp()

    def tearDown(self):
        super(TestEquilibrationSetup, self).tearDown()

    def test_equil(self):
        """ Test for equilibration function."""
        comparetools = testtools.CompareTools(ref_dir, verbose=True)

        if call("python2.7 $PROTOMSHOME/protoms.py -s equilibration -l dcb.pdb -p protein.pdb --nequil 100 --ranseed 100000 --gaff gaff14", shell=True) == 0:

            # Checking whether the required output files have been setup for equilibration MC moves.
            for infile in input_files_setup:
                self.assertTrue(os.path.exists(infile),
                                "Setup input file {0} is missing.".format(infile))

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "Setup output file {0} is missing.".format(outfile))

            print("ProtoMS ligand and protein setup is successful.")

            # Checking content of setup output files with reference data in files.
            for outfile in output_files_setup:
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

        else:
            raise simulationobjects.SetupError("ProtoMS ligand and protein setup failed.")

        if call("$PROTOMSHOME/build/protoms3 run_bnd.cmd", shell=True) == 0:

            # Checking whether the simulation output files have been created successfully for equilibration MC moves.
            for outfile in out_sim_files:
                outfile_rel = os.path.join("out_bnd", outfile)
                self.assertTrue(os.path.exists(outfile_rel),
                                "Simulation file {0} is missing.".format(outfile_rel))

                # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("Simulation was not successful.")

# Entry point to nosetests or unittests.
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
