# This test is to check the equilibration MC moves prior to Sampling MC moves.

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
# ProtoMS Equilibration MC moves test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb", "run_bnd.cmd"]
ref_header_list = ['HEADER', 'cap' , '32.7139', '8.3309', '4.4997', '30.0000', '1.5']
out_sim_files = ["info", "equil_bnd.pdb", "warning"]

class TestEquilibrationSetup(unittest.TestCase):
    
    def setUp(self):
        super(TestEquilibrationSetup, self).setUp()

    def tearDown(self):
        super(TestEquilibrationSetup, self).tearDown()

    def test_equil(self):
        
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s equilibration -l dcb.pdb -p protein.pdb --nequil 100 --ranseed 100000", shell=True)) == 0):

            #Checking whether the required output files have been setup for equilibration MC moves.
                
            for out_files in output_files_setup:
	        try:
                    self.assertTrue(os.path.exists(out_files))
                except IOError as e:
  		    print e
		    print "ProtoMS setup output file ",out_files, "is missing.", "There could be problems with zmat generation, forcefield issues and ProtoMS command files generation."

            #Checking whether the setup output water HEADER is approximately same as reference water.pdb file.

            water_file = open("water.pdb", "r")
 	    header_line = water_file.readline()
            header_line = header_line.replace('\n', '')
            header_list = re.split(" +",header_line)

            #Comparing with reference header
            
            if header_list[0] == ref_header_list[0] and header_list[1] == ref_header_list[1] and \
(header_list[2] == ref_header_list[2]) or ('{:0.2f}'.format(float(header_list[2])) == '{:0.2f}'.format(float(ref_header_list[2]))) and \
(header_list[3] == ref_header_list[3]) or ('{:0.2f}'.format(float(header_list[3])) == '{:0.2f}'.format(float(ref_header_list[3]))) and \
(header_list[4] == ref_header_list[4]) or ('{:0.2f}'.format(float(header_list[4])) == '{:0.2f}'.format(float(ref_header_list[4]))) and \
header_list[5] == ref_header_list[5] and header_list[6] == ref_header_list[6]:
              
                print "ProtoMS ligand and protein setup is successful."

            else:
                print "Discrepancy in HEADER parameters in water cap-file. Please check!" 


	    print "Setup and command files generation is successful."

            if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            #Checking whether the simulation output files have been created successfully for equilibration MC moves
                for out_files in out_sim_files:
                    try:
                        self.assertTrue(os.path.exists("out_bnd/"+ out_files))
                    except IOError as e:
                        print e
                        print "Equilibration simulation file: ",out_sim_files, "is missing."  
                        
            else:
                #logger.error(
                print "Equilibration simulation is not successful."
        else:
            #logger.error(
            print "Equilibration MC move check is is not successful. Either protoms setup or simulation or both failed."
            

if __name__ == '__main__':
    unittest.main()
