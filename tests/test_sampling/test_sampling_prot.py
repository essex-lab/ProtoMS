"""This test is to check the Sampling MC moves."""

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
# ProtoMS Equilibration Sampling MC moves test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/sampling/"

output_files_setup = ["run_bnd.cmd"]
out_sim_files = ["info", "all.pdb", "warning", "results", "restart.prev", "restart", "accept"]


class TestSampling(unittest.TestCase):

    """Test for Sampling function."""

    def setUp(self):
        super(TestSampling, self).setUp()

    def tearDown(self):
        super(TestSampling, self).tearDown()

    def test_sampling(self):
        """ Test for Sampling function."""

        if((call("python2.7 $PROTOMSHOME/protoms.py -s sampling -l dcb.pdb -p protein.pdb --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10", shell=True)) == 0):

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

        if((call("$PROTOMSHOME/build/protoms3 run_bnd.cmd", shell=True)) == 0):
            # Checking whether the simulation output files have been created successfully for Sampling MC moves.
            for outfile in out_sim_files:
                outfile_rel = os.path.join("out_bnd", outfile)
                self.assertTrue(os.path.exists(outfile_rel),
                                "Sampling simulation file {0} is missing.".format(outfile_rel))

                # Checking content of sampling MC moves output files with reference data in files.
                if outfile == "info":
                    if((call("bash content_info_comp.sh", shell=True)) == 0):
                        print "\n Info files content matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference info file.")
                else:
                    self.assertTrue(filecmp.cmp(outfile_rel, os.path.join(ref_dir, outfile_rel)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("Sampling MC protoms simulation failed!")

# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
