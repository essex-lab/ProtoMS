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

class TestFortran(unittest.TestCase):

# Test if ProtoMS Fortran code is compiled correctly.
    def setUp(self):
        super(TestFortran, self).setUp()
        print("Setting up Fortran code compilation test.")
    
    def tearDown(self):
        super(TestFortran, self).setUp()
        print("Cleaning up ProtoMS Fortran code test.")

    def test_build_folder(self):
    
        try:
            os.mkdir(os.environ["PROTOMSHOME"] + "/build")
        except OSError as e:
            print e       
            if ((call("cd $PROTOMSHOME/build && sudo make install", shell = True)) == 0):
                print("Fortran code is compiled correctly and ProtoMS is built.")

    def test_fortran(self):
        

        self.assertTrue(os.path.exists(os.environ["PROTOMSHOME"] + "/protoms3"), "protoms3 executable doesn't exist. Fortran code is not compiled properly.")


#Entry point for nosetests or unittests
	    

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()

