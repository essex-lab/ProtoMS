#!/usr/bin/env python2.7

"""
Tests for checking whether MPI is installed correctly.
"""

import unittest
import nose
import nose.tools as nt
import subprocess

from subprocess import call


class TestMPIInstallation(unittest.TestCase):
    """Check that MPI is installed correctly."""

    def test_mpi_install(self):
        if (call("mpirun --version", shell=True)) != 0 or (call("mpiexec --version", shell=True)) != 0:
            print("No MPI implementation is present.")
            raise ValueError("mpirun or mpiexec is not present")


# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
