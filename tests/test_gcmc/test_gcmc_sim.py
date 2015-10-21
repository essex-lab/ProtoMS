# This test is to check the protoms.py setup and protoms3 run of gcmc simulations.

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
out_gcmc_tools = ["gcmc_box.pdb"]
output_files_setup = ["gcmc_wat.pdb","run_bnd.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning", "info"]

class TestGCMC(unittest.TestCase):
    
    def setUp(self):
        super(TestGCMC, self).setUp()

    def tearDown(self):
        super(TestGCMC, self).tearDown()

    def test_gcmc(self):
        
        if((call("python2.7 $PROTOMSHOME/tools/make_gcmcbox.py -s" + test_dir + "wat.pdb", shell=True)) == 0):
            
            #Checking whether GCMC tool created the box successfully.
            
            for out_files in out_gcmc_tools:
                try:
                    os.path.exists(out_files)
                except IOError as e:
                    print e
                    print "GCMC box file ",test_dir + out_files, "is missing." 
        
            if((call("python2.7 $PROTOMSHOME/protoms.py -sc protein.dcb -s gcmc --gcmcwater wat.pdb --gcmcbox gcmc_box.pdb --adams 20 --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 --capradius 26 -w water.pdb -f" + test_dir, shell=True)) == 0):

                #Checking whether the required output files have been setup for GCMC simulation.
                
                for out_files in output_files_setup:
	            try:
                        os.path.exists(test_dir + out_files)
                    except IOError as e:
  		        print e
		        print "ProtoMS setup output file ",test_dir + out_files, "is missing.", "There could be problems with ProtoMS input command file generation for simulation."


  	        print "Setup and command file generation is successful."

                if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

                #Checking whether the simulation output files have been created successfully for Sampling MC moves
                    for out_files in out_sim_files:
                        try:
                            os.path.exists("out/"+ out_files)
                        except IOError as e:
                            print e
                            print "GCMC simulation file: ",out_files, "is missing."  
                        
                else:
                    print "GCMC simulation is not successful."
        else:
            print "GCMC is not successful. Either protoms setup or simulation or both failed."
            

if __name__ == '__main__':
    unittest.main()
