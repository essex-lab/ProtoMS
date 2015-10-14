# This test is to check the protoms.py Simulation None Setup.

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
# ProtoMS Simulation None Setup test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
test_dir = proto_env + "/tests/test_setup/"
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb"]
ref_header_list = ['HEADER', 'cap' , '32.7139', '8.3309', '4.4997', '30.0000', '1.5']

class TestProtSetup(unittest.TestCase):
    
    def setUp(self):
        super(TestProtSetup, self).setUp()

    def tearDown(self):
        super(TestProtSetup, self).tearDown()

    def test_prep(self):
        
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s none -l dcb.pdb -p protein.pdb -f " + test_dir, shell=True)) == 0):

            #Checking whether the required output files have been setup.
                
            for out_files in output_files_setup:
                try:
                    self.assertTrue(os.path.exists(test_dir + out_files))
                except IOError as e:
                    print e
		    print("ProtoMS setup output file ",output_files_setup, "is missing.", "There could be problems with zmat generation or forcefield issues.")

            #Checking whether the setup output water HEADER is approximately same as reference water.pdb file.

            water_file = open(test_dir+"water.pdb", "r")
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

if __name__ == '__main__':
    unittest.main()
