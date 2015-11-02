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

from subprocess import call


class TestMPIInstallation(unittest.TestCase):

    """Check that MPI is installed correctly."""

    def test_mpi_install(self):

        if (call("mpirun --version", shell=True)) != 0 or (call("mpiexec --version", shell=True)) != 0:
            print "No MPI implementation is present."
            raise ValueError("mpirun or mpiexec is not present")
        else:
            call("mpiexec --version", shell=True)


# Entry point for nosetests or unittests


if  __name__ == '__main__':
    unittest.main()
    nose.runmodule()
          
