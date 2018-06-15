#!/usr/bin/env python

"""
Tests for ProtoMS python module (protoms.py): Setup tests
"""

import os
import unittest
import nose
import nose.tools as nt

# ---------------------------------------------
# PROTOMSHOME test
# ---------------------------------------------


class TestSetUp(unittest.TestCase):
    """Test if PROTOMSHOME is set properly. PROTOMSHOME must be set after protoms install at HOME directory."""

    def setUp(self):
        super(TestSetUp, self).setUp()
        print("Setting up PROTOMSHOME check test.")

    def tearDown(self):
        super(TestSetUp, self).tearDown()
        print("Cleaning up PROTOMSHOME check test.")

    def test_protoms_path(self):
        self.assertIsNotNone(os.getenv("PROTOMSHOME"), "PROTOMSHOME is not set.")


if __name__ == '__main__':
    unittest.main()
