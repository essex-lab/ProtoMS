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
# ProtoMS parameters test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ff_tem_files = ["amber99-residues-mod.ff", "amber99-residues.ff", "amber99.ff", "amber99SB.ff", "gaff.ff", "gaff.types", "gaff14.ff", "gaff14.types", "gborn.parameters", "gbornrevised.parameters","gold-residues.ff", "opls96-residues.ff", "opls96.ff", "rotamlib", "solvents.ff", "surface.parameters"]

class TestParamSetUp(unittest.TestCase):

# Test if ProtoMS reference files and templates exist are in the expected place.
    def setUp(self):
        super(TestToolsSetUp, self).setUp()
        print("Setting up ProtoMS parameters test.")
    
    def tearDown(self):
        super(TestToolsSetUp, self).setUp()
        print("Cleaning up ProtoMS parameters test.")

    def test_params(self):
        
        for tem_files in ff_tem_files:

            try:
                self.assertTrue(os.path.exists(proto_env + "/parameter/" + tem_files))
    
            except IOError as e:
                print e
	        print("Reference files/ force field parameter file ", tem_files, "is not present in the expected place.")


if __name__ == '__main__':
    unittest.main()

