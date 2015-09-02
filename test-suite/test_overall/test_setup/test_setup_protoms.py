#!/usr/bin/python2.7

"""
Tests for ProtoMS python module (protoms.py) setup of protein and ligand.
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
import argparse
import subprocess
import logging
import time

import numpy as np
import re

from protoms import _get_prefix, _locate_file, _merge_templates, _load_ligand_pdb, _prep_ligand, _prep_protein, _prep_singletopology, _prep_gcmc, _prep_jaws2, _cleanup, _wizard

from tools import simulationobjects

#---------------------------------------------
# Setup tests
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$', '', proto_path)


class TestProtLigSetUp(unittest.TestCase):


    def setUp(self):
        """ Initialise test suite - No op. """
        
        super(TestProtLigSetUp, self).setUp()
        print("Setting up protoMS.py Python module test cases.")
        pass

    def tearDown(self):
        """ Clean up test suite - No op. """

        super(TestProtLigSetUp, self).tearDown()
        print("Cleaning up ProtoMS's python module test cases.")
        pass

    def test_locate_file_none(self):
        
        with self.assertRaises(IOError):
            raise simulationobjects.SetUpError("File doesn't exist.")
 
   # def test_protoms(self):
   #     self.assertEqual()                

if  __name__ == '__main__':
    unittest.main()
