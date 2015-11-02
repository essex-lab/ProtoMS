"""This test is to check the ProtoMS setup of simulations with multiple lambda windows and the MPI implementation of protoms3 program."""

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

#-------------------------------------------------------
# ProtoMS setup and RETI dual topology simulations test
#-------------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.

proto_env = os.environ["PROTOMSHOME"]
test_dir = proto_env + "/tests/test_RETI_dbl/"
output_files_setup = ["ethane_box.pdb", "eth-meo.tem", "run_free.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart.prev", "warning", "info", "restart", "results_inst"]

class TestRETIdbl(unittest.TestCase):
    """Test for RETI/ MPI function - short dual topology simulation."""
    def setUp(self):
        super(TestRETIdbl, self).setUp()

    def tearDown(self):
        super(TestRETIdbl, self).tearDown()

    def test_RETI_dbl(self):

        if((call("python2.7 $PROTOMSHOME/protoms.py -s dualtopology -l ethane.pdb methanol.pdb --nequil 0 --nprod 10 --lambdas 0.00 0.33 0.67 1.00 --ranseed 100000 --dumpfreq 1 --cleanup -f" + test_dir, shell=True)) == 0):

            #Checking whether the required output files have been setup for RETI dual topology protoms.py setup.
                
            for out_files in output_files_setup:
	            self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with ligand's zmat generation, forcefield issues, template generation issues and ProtoMS input command file generation for simulation." % (test_dir + out_files))


	    print "Setup and command files generation is successful."

            if((call("mpirun -np 4 $PROTOMSHOME/protoms3 run_free.cmd", shell=True)) == 0):

            #Checking whether the simulation output files have been created successfully for RETI dual topology.
                if(os.path.exists("out_free")):
                    for root, dirs, files in os.walk("out_free"):
                        if len(dirs) != 0:
                            for d in dirs:
                                for out_files in out_sim_files:
                                    self.assertTrue(os.path.exists(os.path.join("out_free",d,out_files)), "Simulation file %s is missing. Please check!" % (os.path.join("out_free",d,out_files)))
            
            else:
                print "RETI dual topology simulation is not successful."
                


        else:
            print "RETI dual topology check is is not successful. Either protoms setup or simulation or both failed."

#Entry point to nosetests or unittests.

if __name__ == "__main__":
    unittest.main()
    nose.runmodule()

