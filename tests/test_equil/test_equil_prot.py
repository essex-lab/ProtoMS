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
test_dir = proto_env + "/tests/test_equil/"
output_files_setup = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "water.pdb", "run_bnd.cmd"]
ref_header_list = ['HEADER', 'cap' , '32.7139', '8.3309', '4.4997', '30.0000', '1.5']
out_sim_files = ["info", "equil_bnd.pdb", "warning"]
outfiles = ["dcb.prepi", "dcb.frcmod", "dcb.zmat", "dcb.tem", "dcb_box.pdb", "protein_scoop.pdb", "run_bnd.cmd"]

class TestEquilibrationSetup(unittest.TestCase):
    
    """Test for equilibration function."""
    
    def setUp(self):
        super(TestEquilibrationSetup, self).setUp()

    def tearDown(self):
        super(TestEquilibrationSetup, self).tearDown()

    def test_equil(self):
        
        """ Test for equilibration function."""
        
        if((call("python2.7 $PROTOMSHOME/protoms.py -s equilibration -l dcb.pdb -p protein.pdb --nequil 100 --ranseed 100000 -f" + test_dir, shell=True)) == 0):

            #Checking whether the required output files have been setup for equilibration MC moves.
                
            for out_files in output_files_setup:
	        
                self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with zmat generation, forcefield issues and ProtoMS command files generation." % (test_dir + out_files))
           
            #Checking whether the setup output water HEADER is approximately same as reference water.pdb file.

            water_file = open(test_dir + "water.pdb", "r")
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

            #Checking content of setup output files with reference data in files.

            for out_files in outfiles:
                if((call("diff "+ test_dir + out_files + " $PROTOMSHOME/tests/equil/" + out_files, shell=True)) == 0):
                    continue
                else:
                    raise ValueError("Content mismatch between output and reference %s" %(out_files))

        else:
            raise simulationobjects.SetupError("ProtoMS ligand and protein setup is not successful.")
    
	   

        if ((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            """Checking whether the simulation output files have been created successfully for equilibration MC moves."""
            for out_files in out_sim_files:
                   
                self.assertTrue(os.path.exists("out_bnd/"+ out_files),"Equilibration simulation file: %s is missing." %("out_bnd/"+ out_files))

            #Checking content of equilibration simulation output files with reference data in files.

                if out_files == "info":
                    if((call("bash content_info_comp.sh", shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference info file.")
                else:
                    if((call("diff "+ test_dir + "out_bnd/" + out_files + " $PROTOMSHOME/tests/equil/out_bnd/" + out_files, shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(os.path.join("$PROTOMSHOME/tests/equil/out_bnd/",out_files)))

            
        else:
            raise simulationobjects.SetupError("Equilibration simulation is not successful.")

#Entry point to nosetests or unittests.

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
