"""This test is to check the protoms.py Simulation None Setup."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# ---------------------------------------------
# ProtoMS Simulation None Setup test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]

site.addsitedir(proto_env)
from tools import simulationobjects

ref_dir = proto_env + "/tests/setup/"
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb"]


class TestProtSetup(unittest.TestCase):
    """ Test for ProtoMS setup function."""

    def setUp(self):
        super(TestProtSetup, self).setUp()

    def tearDown(self):
        super(TestProtSetup, self).tearDown()

    def test_prep(self):
        """ Test for ProtoMS setup function."""

        if call("python2.7 $PROTOMSHOME/protoms.py -s none -l dcb.pdb -p protein.pdb --setupseed 100000", shell=True) == 0:

            # Checking whether the required output files have been setup.
            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "Setup output file {0} is missing.".format(outfile))

            # Checking content of output files with reference data files
            for outfile in output_files_setup:
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference {0}".format(outfile))
        else:
            raise simulationobjects.SetupError("ProtoMS ligand and protein setup is not successful.")


# Entry point to nosetests or unittests.
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
