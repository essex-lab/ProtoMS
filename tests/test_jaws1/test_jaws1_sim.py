"""This test is to check the ProtoMS setup and protoms3 run of JAWS stage 1 simulations."""

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
test_dir = proto_env + "/tests/test_jaws1/"
ref_dir = proto_env + "/tests/jaws1/"

#Storing present working directory path to a variable.
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$','',proto_path)

output_files_setup = ["fragment.tem", "fragment.frcmod", "fragment.prepi", "fragment.zmat", "fragment_box.pdb", "protein_scoop.pdb"]
outfiles_setup = ["water_clr.pdb","jaws1_box.pdb", "jaws1_wat.pdb", "run_bnd.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart", "warning", "info"]


class TestJAWS1(unittest.TestCase):
    
    """Test for JAWS1 function."""
    def setUp(self):
        super(TestJAWS1, self).setUp()

    def tearDown(self):
        super(TestJAWS1, self).tearDown()

    def test_jaws(self):
        
        """Test for JAWS1 function."""
        if((call("python2.7 $PROTOMSHOME/protoms.py -s jaws1 -p protein.pdb -l fragment.pdb --nequil 0 --nprod 100 --ranseed 100000 --dumpfreq 10 -w water.pdb -f" + test_dir, shell=True)) == 0):

            """Checking whether the required output files have been setup for JAWS Stage 1 simulations."""
                
            for out_files in output_files_setup:

	        self.assertTrue(os.path.exists(test_dir + out_files), "ProtoMS setup output file %s is missing. There could be problems with zmat generation, forcefield issues and ProtoMS input command file generation for simulation." % (test_dir + out_files))

            for out_files in outfiles_setup:
                self.assertTrue(os.path.exists(proto_path+"/"+out_files), "ProtoMS setup output file %s is missing. Please check!" % (os.path.join(proto_path,"/",out_files)))


            """ Checking content of setup output files and reference data."""
            for out_files in output_files_setup:
                if out_files == "protein_scoop.pdb":
                    if((call("bash "+test_dir+"content_ps_comp.sh", shell=True)) == 0):
                        print "\n Protein scoop files matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference Protein scoop files.")
                else:
                    if((call("diff " + test_dir + out_files + " " + ref_dir + out_files, shell=True)) == 0):
                        print "\n Content matched for %s." %out_files
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(out_files))

            for out_files in outfiles_setup:
                if out_files == "run_bnd.cmd":
                    if((call("bash "+ test_dir+"content_cmd_comp.sh", shell=True)) == 0):
                        print "\n Output Command files matched."
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference command files.")
                else:
                    if((call("diff "+ out_files + " "+ ref_dir + out_files, shell=True)) == 0):
                        continue
                    else:
                        raise ValueError("Content mismatch between output and reference %s" %(out_files))                  

	else:
            raise simulationobjects.SetupError("ProtoMS setup and command files generation failed!")

        if((call("$PROTOMSHOME/protoms3 run_bnd.cmd", shell=True)) == 0):

            """Checking whether the simulation output files have been created successfully for JAWS Stage 1."""
            for out_files in out_sim_files:
                self.assertTrue(os.path.exists("out/"+ out_files), "JAWS Stage 1 simulation file: %s is missing." % out_files )

                """Comparing content of JAWS Stage 1 simulation output files with reference data."""
                if out_files == "info":
                    if((call("bash "+test_dir+"content_info_comp.sh", shell=True)) == 0):
                        print "\n Info file content matched."
                        continue
                    else:
                        if((call("diff out/"+out_files+" "+ref_dir + "out/"+out_files, shell=True)) == 0):
                            print "\n Content matched for %s." %"out/"+out_files
                            continue
                        else:
                            raise ValueError("Content mismatch between output and reference %s" %(os.path.join("out/", out_files)))
 
                        
        else:
            raise simulationobjects.SetupError("JAWS Stage 1 simulation is not successful.")

#Entry point to unittests or nosetests

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
