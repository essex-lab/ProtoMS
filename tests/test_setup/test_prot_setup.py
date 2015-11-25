"""This test is to check the protoms.py Simulation None Setup."""

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

# Storing present working directory path to a variable.
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$','',proto_path)

test_dir = proto_env + "/tests/test_setup/"
ref_dir = proto_env + "/tests/setup/"
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb"]
ref_header_list = ['HEADER', 'cap' , '32.7139', '8.3309', '4.4997', '30.0000', '1.5']

outfiles = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb"]

class TestProtSetup(unittest.TestCase):
    
    """ Test for ProtoMS setup function."""
    
    def setUp(self):
        super(TestProtSetup, self).setUp()

    def tearDown(self):
        super(TestProtSetup, self).tearDown()

    def test_prep(self):
        
        """ Test for ProtoMS setup function."""
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s none -l dcb.pdb -p protein.pdb -f" + test_dir, shell=True)) == 0):

            #Checking whether the required output files have been setup.
                
            for out_files in output_files_setup:
                
                if out_files == "water.pdb":
                    self.assertTrue(os.path.exists(proto_path +"/"+"water.pdb"), "Water-cap file is missing.")
                else:
                    self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file (%s) is missing. There could be problems with zmat generation or forcefield issues." % out_files)
            
            #Checking whether the setup output water HEADER is approximately same as reference water.pdb file.

            water_file = open(proto_path+"/"+"water.pdb", "r")
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
                raise ValueError("Discrepancy in HEADER parameters in water cap-file. Please check!")

            #Checking content of output files with reference data files

            for out_files in outfiles:

                if out_files == "protein_scoop.pdb":
                    if((call("bash "+test_dir+"content_ps_comp.sh", shell=True)) == 0):
                        print "\n Protein scoop files matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference Protein scoop files.")
                else:
                    if((call("diff "+test_dir + out_files + " " +ref_dir + out_files, shell=True)) == 0):
                        print "\n Contents matched for %s." %out_files
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(out_files))
        else:
            raise simulationobjects.SetupError("ProtoMS ligand and protein setup is not successful.")


#Entry point to nosetests or unittests.

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
