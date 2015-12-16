"""This test is to check the protoms.py Simulation None Setup."""

import nose
import unittest
import argparse
import os
import sys
import subprocess
import logging
import time
import re
import numpy as np
import filecmp
import protoms

from protoms import _is_float, _get_prefix, _locate_file, _merge_templates
from protoms import _load_ligand_pdb, _prep_ligand, _prep_protein, _prep_singletopology
from protoms import _prep_gcmc, _prep_jaws2, _cleanup, _wizard

import tools
from tools import simulationobjects

from subprocess import call

# ---------------------------------------------
# ProtoMS Simulation None Setup test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]

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

        if((call("python2.7 $PROTOMSHOME/protoms.py -s none -l dcb.pdb -p protein.pdb --setupseed 100000", shell=True)) == 0):

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
