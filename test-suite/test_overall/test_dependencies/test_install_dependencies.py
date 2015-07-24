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

amber_env = os.environ["AMBERHOME"]
amber_path = amber_env + "/AmberTools/bin"

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

#Checking for AmberTools installation dependencies

     def test_ambertools_dependencies(self):
      
         try:
             self.assertTrue(os.getenv("AMBERHOME"))
         except KeyError as e:
             print e
    	     print("AMBERHOME is not set.")

         try:
             self.assertTrue(os.path.exists(amber_path + "/antechamber"))
	 except IOError as e:
  	     print e
  	     print("antechamber AmberTools module doesn't exist.")
    
	 try:
             self.assertTrue(os.path.exists(amber_path + "/parmchk"))
	 except IOError as e:
  	     print e
  	     print("parmchk AmberTools module doesn't exist.")


if  __name__ == '__main__':
    unittest.main()

          
