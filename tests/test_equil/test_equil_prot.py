"""This test is to check the equilibration MC moves prior to Sampling MC moves."""

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
# ProtoMS Equilibration MC moves test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]

# Storing present working directory path to a variable.
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$', '', proto_path)

test_dir = proto_env + "/tests/test_equil/"
ref_dir = proto_env + "/tests/equil/"
input_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem",
                     "dcb_box.pdb", "protein_scoop.pdb", "water.pdb"]
output_files_setup = ["run_bnd.cmd"]
out_sim_files = ["info", "equil_bnd.pdb", "warning"]


class TestEquilibrationSetup(unittest.TestCase):
    """Test for equilibration function."""

    def setUp(self):
        super(TestEquilibrationSetup, self).setUp()

    def tearDown(self):
        super(TestEquilibrationSetup, self).tearDown()

    def test_equil(self):
        """ Test for equilibration function."""

        if((call("python2.7 $PROTOMSHOME/protoms.py -s equilibration -l dcb.pdb -p protein.pdb --nequil 100 --ranseed 100000", shell=True)) == 0):

            # Checking whether the required output files have been setup for equilibration MC moves.
            for infile in input_files_setup:
                self.assertTrue(os.path.exists(test_dir + infile),
                                "Setup input file {0} is missing.".format(infile))

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(test_dir + outfile),
                                "Setup output file {0} is missing.".format(outfile))

            print "ProtoMS ligand and protein setup is successful."

            # Checking content of setup output files with reference data in files.
            for outfile in output_files_setup:
                self.assertTrue(filecmp.cmp(test_dir+outfile, ref_dir+outfile),
                                "Content mismatch between output and reference for file {0}".format(outfile))

        else:
            raise simulationobjects.SetupError("ProtoMS ligand and protein setup failed.")

        if ((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            # Checking whether the simulation output files have been created successfully for equilibration MC moves.
            for outfile in out_sim_files:
                outfile_rel = os.path.join("out_bnd", outfile)
                self.assertTrue(os.path.exists(outfile_rel),
                                "Simulation output file {0} is missing.".format(outfile_rel))

            # Checking content of equilibration simulation output files with reference data in files.
                if outfile == "info":
                    if((call("bash "+test_dir+"content_info_comp.sh", shell=True)) == 0):
                        print "\n Info file contents matched."
                    else:
                        raise ValueError("Content mismatch between output and reference info file.")

                else:
                    self.assertTrue(filecmp.cmp(outfile_rel, os.path.join(ref_dir, outfile_rel)),
                                    "Content mismatch between output and reference for file {0}".format(os.path.join("out_bnd", outfile)))

        else:
            raise simulationobjects.SetupError("Simulation was not successful.")

# Entry point to nosetests or unittests.
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
