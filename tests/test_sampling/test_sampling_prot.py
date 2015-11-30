"""This test is to check the Sampling MC moves."""

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
# ProtoMS Equilibration Sampling MC moves test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
test_dir = proto_env + "/tests/test_sampling/"
ref_dir = proto_env + "/tests/sampling/"

#Storing present working directory path to a variable.
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$','',proto_path)

ref_header_list = ['HEADER', 'cap' , '32.7139', '8.3309', '4.4997', '30.0000', '1.5']
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb", "run_bnd.cmd"]
out_sim_files = ["info", "all.pdb", "warning", "results", "restart.prev", "restart", "accept"]
outfiles = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "run_bnd.cmd"]

class TestSampling(unittest.TestCase):
    
    """Test for Sampling function."""

    def setUp(self):
        super(TestSampling, self).setUp()

    def tearDown(self):
        super(TestSampling, self).tearDown()

    def test_sampling(self):
        
        """ Test for Sampling function."""
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s sampling -l dcb.pdb -p protein.pdb --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 -f" + test_dir, shell=True)) == 0):

            """Checking whether the required output files have been setup for Sampling MC moves."""
                
            for out_files in output_files_setup:
                
                if out_files == "water.pdb" or out_files == "run_bnd.cmd":
                    self.assertTrue(os.path.exists(proto_path +"/"+out_files), "%s file is missing." %out_files)
                else:

                    self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with zmat generation, forcefield issues and ProtoMS input command file generation for simulation." % (test_dir + out_files))

                
            for out_files in output_files_setup:

                """Checking content of setup output files with reference data in files."""
                if out_files == "protein_scoop.pdb":
                    if((call("bash "+test_dir+"content_ps_comp.sh", shell=True)) == 0):
                        print "\n Protein scoop files matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference Protein scoop files.")

                elif out_files == "water.pdb":

                    #Checking whether the setup output water HEADER is approximately same as reference water.pdb file.

                    water_file = open(proto_path + "/" + "water.pdb", "r")
 	            header_line = water_file.readline()
                    header_line = header_line.replace('\n', '')
                    header_list = re.split(" +",header_line)

                    #Comparing with reference header
            
                    if header_list[0] == ref_header_list[0] and header_list[1] == ref_header_list[1] and \
(header_list[2] == ref_header_list[2]) or ('{:0.2f}'.format(float(header_list[2])) == '{:0.2f}'.format(float(ref_header_list[2]))) and \
(header_list[3] == ref_header_list[3]) or ('{:0.2f}'.format(float(header_list[3])) == '{:0.2f}'.format(float(ref_header_list[3]))) and \
(header_list[4] == ref_header_list[4]) or ('{:0.2f}'.format(float(header_list[4])) == '{:0.2f}'.format(float(ref_header_list[4]))) and \
header_list[5] == ref_header_list[5] and header_list[6] == ref_header_list[6]:
              
                        print("Water-cap file headers matched.")
                        continue
                    else:
                        raise ValueError("Discrepancy in HEADER parameters in water cap-file. Please check!")

                
                elif out_files == "run_bnd.cmd":
                    if((call("bash "+ test_dir+"content_cmd_comp.sh", shell=True)) == 0):
                        print "\n Output Command files matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference command files.")


                else:
                    if((call("diff "+ test_dir + out_files + " "+ ref_dir + out_files, shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(out_files))
  
        else:
            raise simulationobjects.SetupError("ProtoMS setup for sampling MC moves is not successful!")

        if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            """Checking whether the simulation output files have been created successfully for Sampling MC moves."""
            for out_files in out_sim_files:
                self.assertTrue(os.path.exists("out_bnd/"+ out_files), "Sampling simulation file: %s is missing." % out_files )

                """Checking content of sampling MC moves output files with reference data in files."""
        
                if out_files == "info":
                    if((call("bash "+ test_dir+"content_info_comp.sh", shell=True)) == 0):
                        print "\n Info files content matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference info file.")
                else:
                    if((call("diff out_bnd/" + out_files + " "+ ref_dir+"out_bnd/" + out_files, shell=True)) == 0):
                        print "\n %s content matched." %"out_bnd/" + out_files
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(os.path.join("out_bnd/", out_files)))
            
        else:
            raise simulationobjects.SetupError("Sampling MC protoms simulation failed!")

#Entry point for nosetests or unittests

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()

