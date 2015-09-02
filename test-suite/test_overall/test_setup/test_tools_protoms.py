#!/usr/bin/python2.7

"""
Tests for ProtoMS python module (protoms.py): Setup tests
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
from protoms import *

#---------------------------------------------
# Setup tests
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
proto_path = os.environ["HOME"] + "/protoms"

class TestToolsSetUp(unittest.TestCase):

# Test if ProtoMS tools and reference files exist are in the expected place.
    def setUp(self):
        super(TestToolsSetUp, self).setUp()
        print("Setting up ProtoMS tools test.")
    
    def tearDown(self):
        super(TestToolsSetUp, self).setUp()
        print("Cleaning up ProtoMS tools test.")

    def setup_test_tools(self):

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/ambertools.py"))        

        except:
            print("AmberTools doesn't exist. antechamber and parmchk can't be encapsulated from AmberTools suite of programs.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/build_template.py"))

        except:
            print("ProtoMS template build tool doesn't exist.")

	
 	try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_bar.py"))

	except:
	    print("Free energy calculation tool doesn't exist.")
             
 
        try:
	    self.assertTrue(os.path.isfile(proto_path + "/tools/calc_clusters.py"))

        except:
            print("Molecular clustering tool doesn't exist.")
    

	try:
	    self.assertTrue(os.path.isfile(proto_path + "/tools/calc_density.py"))

        except:
            print("Density calculation tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_dg.py"))

        except:
            print("Free energy calculation (using TI, BAR and MBAR) tool doesn't exist.")


	try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_gcsingle.py"))

        except:
            print("Free energy calculation and analysis tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_replicapath.py"))

        except:
            print("Replica plotting tool doesn't exist.")	

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_rmsd.py"))

        except:
            print("Ligand RMSD calculation tool doesn't exist.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_series.py"))

        except:
            print("Time series plotting and analysing tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/calc_ti.py"))

        except:
            print("Free energy calculation (TI) tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.isfile(proto_path +	"/tools/clear_gcmcbox.py"))

        except:
            print("Simulation box clearing tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/convertatomnames.py"))

        except:
            print("Naming conversion tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/convertwater.py"))

        except:
            print("Water molecule conversion tool doesn't exist.")

 
        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/distribute_waters.py"))

        except:
            print("Water molecules placing tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/divide_pdb.py"))

        except:
            print("PDB file splitting tool doesn't exist.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/generate_input.py"))

        except:
            print("Input files generating tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/make_dummy.py"))

        except:
            print("Matching dummy particle maker tool doesn't exist.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/make_gcmcbox.py"))

        except:
            print("GCMC/ JAWS simulation box maker tool doesn't exist.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/make_single.py"))

        except:
            print("ProtoMS template files generator for single topology free energy simulations doesn't exist")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/merge_templates.py"))

        except:
            print("ProtoMS template combining tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/plot_theta.py"))

        except:
            print("Theta distribution plotting tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/pms2pymbar.py"))
        except:
            print("Free energy extraction tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/scoop.py"))
        except:
            print("Protein truncating tool doesn't exist.")


        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/solvate.py"))
        except:
            print("Ligand solvating tool doesn't exist.")

        try:
            self.assertTrue(os.path.isfile(proto_path + "/tools/split_jawswater.py"))
        except:
            print("PDB file splitting tool doesn't exist.")   


if  __name__ == '__main__':
    unittest.main()
