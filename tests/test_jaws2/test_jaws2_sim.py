"""This test is to check the ProtoMS setup and protoms3 run of JAWS stage 2 simulations."""

import nose
import unittest
import argparse
import os
import sys
import subprocess
import logging
import time
import numpy as np
import protoms

from protoms import _is_float, _get_prefix, _locate_file, _merge_templates, _load_ligand_pdb, _prep_ligand, _prep_protein, _prep_singletopology, _prep_gcmc, _prep_jaws2, _cleanup, _wizard

import tools
from tools import simulationobjects

from subprocess import call

#------------------------------------------------
# ProtoMS setup and JAWS stage 2 simulations test
#------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.

proto_env = os.environ["PROTOMSHOME"]
test_dir = proto_env + "/tests/test_jaws2/"
output_files_setup = ["fragment.tem", "fragment.frcmod", "fragment.prepi", "fragment.zmat", "fragment_box.pdb", "protein_scoop.pdb", "jaws2_wat1.pdb", "jaws2_wat2.pdb", "jaws2_wat3.pdb", "jaws2_not1.pdb", "jaws2_not2.pdb", "jaws2_not3.pdb", "run_jaws2-w1_bnd.cmd", "run_jaws2-w2_bnd.cmd", "run_jaws2-w3_bnd.cmd"]
out_sim_files = ["accept", "all.pdb", "info", "restart", "restart.prev", "results", "warning"]

class TestJAWS2(unittest.TestCase):
    
    """Test for JAWS2 function."""
    def setUp(self):
        super(TestJAWS2, self).setUp()

    def tearDown(self):
        super(TestJAWS2, self).tearDown()

    def test_jaws2(self):
        
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s jaws2 -l fragment.pdb -p protein.pdb --gcmcwater " + test_dir + "jaws2_waters.pdb --jawsbias 8 10 12 14 --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 -f" + test_dir, shell=True)) == 0):

            #Checking whether the required output files have been setup for JAWS Stage 2 protoms.py setup.
                
            for out_files in output_files_setup:
	            self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with zmat generation, forcefield issues and ProtoMS input command file generation for simulation." % (test_dir + out_files))

	    print "Setup and command files generation is successful."

            if((call("mpirun -np 4 $PROTOMSHOME/protoms3 run_jaws2-w1_bnd.cmd", shell=True)) == 0):

            #Checking whether the simulation output files have been created successfully for JAWS Stage 2.
                if(os.path.exists("out_jaws2-w1")):
                    for root, dirs, files in os.walk("out_jaws2-w1"):
                        if len(dirs) != 0:
                            for d in dirs:
                                for out_files in out_sim_files:
                                    self.assertTrue(os.path.exists(os.path.join("out_jaws2-w1",d,out_files)), "Simulation file %s is missing. Please check!" % os.path.join("out_jaws2-w1",d,out_files))
            
            else:
                print "JAWS Stage 2 simulation is not successful."
        else:
            print "JAWS Stage 2 check is is not successful. Either protoms setup or simulation or both failed."

#Entry point to unittests or nosetests

if __name__ == "__main__":
    unittest.main()
    nose.runmodule()


