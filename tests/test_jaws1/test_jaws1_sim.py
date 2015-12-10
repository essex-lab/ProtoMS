"""This test is to check the ProtoMS setup and protoms3 run of JAWS stage 1 simulations."""

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

# ------------------------------------------------
# ProtoMS setup and JAWS stage 1 simulations test
# ------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/jaws1/"

output_files_setup = ["fragment.tem", "fragment.frcmod", "fragment.prepi", "fragment.zmat", "fragment_box.pdb", "protein_scoop.pdb"]
outfiles_setup = ["water_clr.pdb","jaws1_box.pdb", "jaws1_wat.pdb", "run_jaws.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning", "info"]


class TestJAWS1(unittest.TestCase):
    """Test for JAWS1 function."""

    def setUp(self):
        super(TestJAWS1, self).setUp()

    def tearDown(self):
        super(TestJAWS1, self).tearDown()

    def test_jaws(self):
        """Test for JAWS1 function."""

        if((call("python2.7 $PROTOMSHOME/protoms.py -s jaws1 -p protein.pdb -l fragment.pdb --nequil 0 --nprod 100 --ranseed 100000 --setupseed 100000 --dumpfreq 10 -w water.pdb", shell=True)) == 0):
            # Checking whether the required output files have been setup for JAWS Stage 1 simulations.

            for outfile in output_files_setup:
	        self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            for outfile in outfiles_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            # Checking content of setup output files and reference data.
            for outfile in output_files_setup:
                if outfile == "protein_scoop.pdb":
                    self.assertTrue(call("bash content_ps_comp.sh", shell=True) == 0,
                                    "Content mismatch between output and reference for file {0}".format(outfile))
                else:
                    self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile))

            for outfile in outfiles_setup:
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

	else:
            raise simulationobjects.SetupError("ProtoMS setup and command files generation failed!")

        if((call("$PROTOMSHOME/build/protoms3 run_jaws.cmd", shell=True)) == 0):

            # Checking whether the simulation output files have been created successfully for JAWS Stage 1.
            for outfile in out_sim_files:
                outfile_rel = os.path.join("out_jaws", outfile)
                self.assertTrue(os.path.exists(outfile_rel),
                                "JAWS Stage 1 simulation file {0} is missing.".format(outfile_rel))

                # Comparing content of JAWS Stage 1 simulation output files with reference data.
                if outfile == "info":
                    self.assertTrue(call("bash content_info_comp.sh", shell=True) == 0,
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))
                else:
                    self.assertTrue(filecmp.cmp(outfile_rel, os.path.join(ref_dir, outfile_rel)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("JAWS Stage 1 simulation was not successful.")

# Entry point to unittests or nosetests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
