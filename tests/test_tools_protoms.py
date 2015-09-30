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
import re

#---------------------------------------------
# Setup tests
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
#proto_path = os.popen("pwd").read() 
#proto_path = re.sub('\\n$', '', proto_path)


class TestToolsSetUp(unittest.TestCase):

# Test if ProtoMS tools and reference files exist are in the expected place.
    def setUp(self):
        super(TestToolsSetUp, self).setUp()
        print("Setting up ProtoMS tools test.")
    
    def tearDown(self):
        super(TestToolsSetUp, self).setUp()
        print("Cleaning up ProtoMS tools test.")

    def test_tools(self):

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/ambertools.py"))        

        except IOError as e:
            print e
            print("AmberTools doesn't exist. antechamber and parmchk can't be encapsulated from AmberTools suite of programs.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/build_template.py"))

        except IOError as e:
            print e
            print("ProtoMS template build tool doesn't exist.")

	
 	try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_bar.py"))

	except IOError as e:
		print e
	    	print("Free energy calculation tool doesn't exist.")
             
 
        try:
	    self.assertTrue(os.path.exists(proto_env + "/tools/calc_clusters.py"))

        except IOError as e:
            print e
            print("Molecular clustering tool doesn't exist.")
    

	try:
	    self.assertTrue(os.path.exists(proto_env + "/tools/calc_density.py"))

        except IOError as e:
            print e
            print("Density calculation tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_dg.py"))

        except IOError as e:
            print e
            print("Free energy calculation (using TI, BAR and MBAR) tool doesn't exist.")


	try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_gcsingle.py"))

        except IOError as e:
            print e
            print("Free energy calculation and analysis tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env+ "/tools/calc_replicapath.py"))

        except IOError as e:
            print e
            print("Replica plotting tool doesn't exist.")	

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_rmsd.py"))

        except IOError as e:
            print e
            print("Ligand RMSD calculation tool doesn't exist.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_series.py"))

        except IOError as e:
            print e
            print("Time series plotting and analysing tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/calc_ti.py"))

        except IOError as e:
            print e
            print("Free energy calculation (TI) tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.exists(proto_env +	"/tools/clear_gcmcbox.py"))

        except IOError as e:
            print e
            print("Simulation box clearing tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/convertatomnames.py"))

        except IOError as e:
            print e
            print("Naming conversion tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/convertwater.py"))

        except IOError as e:
            print e
            print("Water molecule conversion tool doesn't exist.")

 
        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/distribute_waters.py"))

        except IOError as e:	
            print e
            print("Water molecules placing tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/divide_pdb.py"))

        except IOError as e:
            print e
            print("PDB file splitting tool doesn't exist.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/generate_input.py"))

        except IOError as e:
            print e
            print("Input files generating tool doesn't exist.")

        
        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/make_dummy.py"))

        except IOError as e:
            print e
            print("Matching dummy particle maker tool doesn't exist.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/make_gcmcbox.py"))

        except IOError as e:
            print e
            print("GCMC/ JAWS simulation box maker tool doesn't exist.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/make_single.py"))

        except IOError as e:
            print e
            print("ProtoMS template files generator for single topology free energy simulations doesn't exist")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/merge_templates.py"))

        except IOError as e:
            print e
            print("ProtoMS template combining tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/plot_theta.py"))

        except IOError as e:
            print e
            print("Theta distribution plotting tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/pms2pymbar.py"))
        except IOError as e:
            print e
            print("Free energy extraction tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/scoop.py"))
        except IOError as e:
            print e
            print("Protein truncating tool doesn't exist.")


        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/solvate.py"))
        except IOError as e:
            print e
            print("Ligand solvating tool doesn't exist.")

        try:
            self.assertTrue(os.path.exists(proto_env + "/tools/split_jawswater.py"))
        except IOError as e:
            print e
            print("PDB file splitting tool doesn't exist.")   


if  __name__ == '__main__':
    unittest.main()
