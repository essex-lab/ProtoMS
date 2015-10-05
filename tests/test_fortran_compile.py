#!/usr/bin/python2.7

"""
Tests for ProtoMS compiled Fortran routines  (protoms3): Simulation test
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
import subprocess

from subprocess import call
#---------------------------------------------
# ProtoMS simulation test
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
proto_build = proto_env + "/build"

class TestFortran(unittest.TestCase):

# Test if ProtoMS Fortran code is compiled correctly.
    def setUp(self):
        super(TestFortran, self).setUp()
        print("Setting up Fortran code compilation test.")
    
    def tearDown(self):
        super(TestFortran, self).setUp()
        print("Cleaning up ProtoMS Fortran code test.")

    def test_fortran(self):
        
        try:
            self.assertTrue(os.path.exists(proto_env + "/protoms3"))
    
        except IOError as e:
            print e
	    print("protoms3 executable doesn't exist. Fortran code is not compiled properly.")


    def test_build_folder(self):
    
        try:
            os.mkdir(proto_env + "/build")
        except OSError as e:
            print e       
            if ((call("cd $PROTOMSHOME/build && make install", shell = True)) == 0):
                print("Fortran code is compiled correctly and ProtoMS is built.")
	    

if __name__ == '__main__':
    unittest.main()

