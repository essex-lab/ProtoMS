#!/usr/bin/env python2.7

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

    """Test if ProtoMS Fortran code is compiled correctly."""
    
    def setUp(self):
        super(TestFortran, self).setUp()
        print("Setting up Fortran code compilation test.")
    
    def tearDown(self):
        super(TestFortran, self).setUp()
        print("Cleaning up ProtoMS Fortran code test.")

    def test_build_folder(self):
        """Test if ProtoMS Fortran code is compiled correctly."""
    
        if not os.path.exists((os.environ["PROTOMSHOME"] + "/build")):
            os.makedirs(os.environ["PROTOMSHOME"] + "/build")

        if (call("cd $PROTOMSHOME/build && cmake .. && make install || sudo make install", shell = True) == 0):
            print("Fortran code is compiled correctly and ProtoMS is successfully built.")
        else:
            raise Exception("Either build files have not been written successfully or make install failed.")
    

    def test_fortran(self):
        """Test to check if protoms3 executable exists.""" 
        
        self.assertTrue(os.path.exists(os.environ["PROTOMSHOME"] + "/protoms3"), "protoms3 executable doesn't exist. Fortran code is not compiled properly.")


#Entry point for nosetests or unittests
	    

if __name__ == '__main__':
    unittest.main()
    nose.runmodule()

