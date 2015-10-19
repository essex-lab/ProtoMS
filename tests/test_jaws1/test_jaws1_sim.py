# This test is to check the ProtoMS setup and protoms3 run of JAWS stage 1 simulations.

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

#------------------------------------------------
# ProtoMS setup and JAWS stage 1 simulations test
#------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
output_files_setup = ["fragment.tem", "fragment.frcmod", "fragment.prepi", "fragment.zmat", "fragment_box.pdb", "jaws1_wat.pdb", "jaws1_box.pdb", "water_clr.pdb", "run_bnd.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning", "info"]

class TestJAWS1(unittest.TestCase):
    
    def setUp(self):
        super(TestJAWS1, self).setUp()

    def tearDown(self):
        super(TestJAWS1, self).tearDown()

    def test_jaws(self):
        
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s jaws1 -p protein.pdb -l fragment.pdb --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 -w water.pdb", shell=True)) == 0):

            #Checking whether the required output files have been setup for JAWS Stage 1 simulations.
                
            for out_files in output_files_setup:
	        try:
                    self.assertTrue(os.path.exists(out_files))
                except IOError as e:
  		    print e
		    print "ProtoMS setup output file ",out_files, "is missing.", "There could be problems with zmat generation, forcefield issues and ProtoMS input command file generation for simulation."


	    print "Setup and command files generation is successful."

            if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            #Checking whether the simulation output files have been created successfully for JAWS Stage 1.
                for out_files in out_sim_files:
                    try:
                        self.assertTrue(os.path.exists("out/"+ out_files))
                    except IOError as e:
                        print e
                        print "JAWS Stage 1 simulation file: ",out_files, "is missing."  
                        
            else:
                print "JAWS Stage 1 simulation is not successful."
        else:
            print "JAWS Stage 1 check is is not successful. Either protoms setup or simulation or both failed."
            

if __name__ == '__main__':
    unittest.main()
