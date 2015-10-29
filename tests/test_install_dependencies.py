#!/usr/bin/python2.7

"""
Tests for checking whether Python and Ambertools dependencies are installed.
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt


class TestDependInstallation(unittest.TestCase):

#Check that Python dependencies are installed correctly.

     def test_python_dependencies(self):

         import numpy
	
         import scipy
  
         import matplotlib


# Checking for AmberTools installation dependencies

     def test_ambertools_dependencies(self):
      
         self.assertTrue(os.getenv("AMBERHOME"))

         self.assertTrue(os.path.exists(os.environ["AMBERHOME"]+ "/AmberTools/bin/antechamber"))

         self.assertTrue(os.path.exists(os.environ["AMBERHOME"]+ "/AmberTools/bin/parmchk"))


#Entry point for nosetests or unittests

if  __name__ == '__main__':
    unittest.main()
    nose.runmodule()
          
