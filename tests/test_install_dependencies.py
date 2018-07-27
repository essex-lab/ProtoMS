#!/usr/bin/env

"""
Tests for checking whether Python and Ambertools dependencies are installed.
"""

import os
import unittest
import nose


class TestDependInstallation(unittest.TestCase):
    """Check that Python dependencies are installed correctly."""

    def test_python_dependencies(self):
        import numpy
        import scipy
        import matplotlib
        import six
        import pymbar

    def test_ambertools_dependencies(self):
        """Checking for AmberTools installation dependencies."""
        self.assertTrue(os.getenv("AMBERHOME"), "AMBERHOME is not set.")

        self.assertTrue(os.path.exists(os.path.join(os.environ["AMBERHOME"],
                                                    "bin/antechamber")),
                        "Antechamber AMBERTOOLS module doesn't exist.")

        self.assertTrue(os.path.exists(os.path.join(os.environ["AMBERHOME"],
                                                    "bin/parmchk2")),
                        "Parmchk AMBERTOOLS module doesn't exist.")


# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
