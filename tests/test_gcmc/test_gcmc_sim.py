"""This test is to check the protoms.py setup and protoms3 run of gcmc simulations."""

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
import protoms

from protoms import _is_float, _get_prefix, _locate_file, _merge_templates, _load_ligand_pdb, _prep_ligand, _prep_protein, _prep_singletopology, _prep_gcmc, _prep_jaws2, _cleanup, _wizard

import tools
from tools import simulationobjects

from subprocess import call

#---------------------------------------------
# ProtoMS setup and GCMC simulations test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
test_dir = proto_env + "/tests/test_gcmc/"
ref_dir = proto_env + "tests/gcmc/"
out_gcmc_tools = ["gcmc_box.pdb"]
output_files_setup = ["gcmc_wat.pdb","run_bnd.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning", "info"]
outfiles = ["gcmc_box.pdb", "gcmc_wat.pdb","run_bnd.cmd"]

class TestGCMC(unittest.TestCase):
    
    """Test for GCMC function."""
    
    def setUp(self):
        super(TestGCMC, self).setUp()

    def tearDown(self):
        super(TestGCMC, self).tearDown()

    def test_gcmc(self):
        
        """Test for GCMC function."""
        
        if((call("python2.7 $PROTOMSHOME/tools/make_gcmcbox.py -s" + test_dir + "wat.pdb", shell=True)) == 0):
            
            """Checking whether GCMC tool created the box successfully."""
            
            for out_files in out_gcmc_tools:
                self.assertTrue(os.path.exists(out_files), "GCMC box file %s is missing." % (test_dir + out_files))
        else:
            raise simulationobjects.SetupError("Making of Simulation box for GCMC is not successful!")
            
        if((call("python2.7 $PROTOMSHOME/protoms.py -sc protein.dcb -s gcmc --gcmcwater wat.pdb --gcmcbox gcmc_box.pdb --adams 20 --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 --capradius 26 -w water.pdb -f" + test_dir, shell=True)) == 0):

            """Checking whether the required output files have been setup for GCMC simulation."""
                
            for out_files in output_files_setup:
                self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with ProtoMS input command file generation for simulation." % (test_dir + out_files))

  	        print "Setup and command files generation is successful."

            """ Checking content of GCMC box and setup output files with reference data in files. """
            for out_files in outfiles:
                if((call("diff "+ test_dir + out_files + " $PROTOMSHOME/tests/gcmc/" + out_files, shell=True)) == 0):
                    continue
                else:
                    raise ValueError("Content mismatch between output and reference %s" %(out_files))

        else:
            raise simulationobjects.SetupError("ProtoMS setup and command files generation is not successful.")

        if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            """Checking whether the simulation output files have been created successfully."""
            for out_files in out_sim_files:
                self.assertTrue(os.path.exists("out/"+ out_files),"GCMC simulation file: %s is missing." % out_files)

                """ Checking content of GCMC simulation output files with reference data in files."""
                if out_files == "info":
                    if((call("bash content_info_comp.sh", shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference info files.")
                else:
                    if((call("diff "+test_dir+ "out/" + out_files+ " $PROTOMSHOME/tests/gcmc/out/" + out_files, shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(os.path.join(ref_dir, "out/", out_files)))
                        
        else:
            raise simulationobjects.SetupError("GCMC simulation is not successful.")
        

#Entry point for nosetests or unittests

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
