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
         
         try:
             import numpy
	 except ImportError as e:
             print e

         
         try:
             import scipy
         except ImportError as e:
             print e

         try:
             import matplotlib
         except ImportError as e:
             print e

     def test_ambertools_dependencies(self):
      
         try:
             self.assertTrue(os.getenv("AMBERHOME"))
         except KeyError as e:
             print e

         

if  __name__ == '__main__':
    unittest.main()

          
