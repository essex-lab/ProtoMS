#!/usr/bin/env python

"""
Tests for checking whether MPI is installed correctly.
"""

#
# Usage examples:
# 
# $ nosetests test_mpi_install.py 
# OR
# $ python test_mpi_install.py
#
#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
import subprocess


class TestMPIInstallation(unittest.TestCase):

# Check that MPI is installed correctly.


    def test_mpi_install(self):

        if os.popen("mpirun --version").read() == '' or os.popen("mpiexec --version").read() == '':
            print "No MPI implementation is present."
            raise ValueError
        else:
            print os.popen("mpiexec --version").read()

# Entry point for nosetests or unittests


if  __name__ == '__main__':
    unittest.main()
    nose.runmodule()
          
