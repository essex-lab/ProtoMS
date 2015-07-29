#!/usr/bin/python2.7

"""
Tests for checking whether MPI is installed correctly.
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
#import subprocess
#import Execute

mpi_path = "/usr/local/bin"

return_check_res = os.popen("mpicxx -v").read()
return_openmpi_info = os.popen(mpi_path+"/ompi_info").read()

class TestMPIInstallation(unittest.TestCase):


     def test_mpi_install(self):

         if not return_check_res or not return_openmpi_info:
             print("OpenMPI support is not installed.")
         else:              
             # Checking whether mpicc and mpirun are present if openMPI is installed.

             try:
                 self.assertTrue(os.path.exists(mpi_path+"/mpicc"))
             except IOError as e:
                 print e
    	         print("mpicc is not present. Check if openMPI is installed correctly.")

             try:
                 self.assertTrue(os.path.exists(mpi_path+"/mpirun"))
 	     except:
                 print e
                 print("mpirun is not present. Check whether OpenMPI is installed correctly.")
        

if  __name__ == '__main__':
    unittest.main()

          
