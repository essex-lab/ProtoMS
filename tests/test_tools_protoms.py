#!/usr/bin/python2.7

"""
Tests for ProtoMS Tools (protoms.py): Setup tests
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
import re

#---------------------------------------------
# ProtoMS tools tests
#---------------------------------------------

# List containing ProtoMS tools.

protoms_tools = ["ambertools.py", "build_template.py", "calc_bar.py", "calc_clusters.py", "calc_density.py", "calc_dg.py", "calc_gcsingle.py", "calc_replicapath.py", "calc_rmsd.py", "calc_series.py", "calc_ti.py", "clear_gcmcbox.py", "convertatomnames.py", "convertwater.py", "distribute_waters.py", "divide_pdb.py", "generate_input.py", "make_dummy.py", "make_gcmcbox.py", "make_single.py", "merge_templates.py", "plot_theta.py", "pms2pymbar.py", "scoop.py", "solvate.py", "split_jawswater.py" ]

class TestToolsSetUp(unittest.TestCase):

    """Test if ProtoMS tools and reference files exist are in the expected place."""     

    def setUp(self):
        super(TestToolsSetUp, self).setUp()
        print("Setting up ProtoMS tools test.")
    
    def tearDown(self):
        super(TestToolsSetUp, self).setUp()
        print("Cleaning up ProtoMS tools test.")

    def test_tools(self):

        for tools_files in protoms_tools:
            self.assertTrue(os.path.exists(os.environ["PROTOMSHOME"] + "/tools/"+ tools_files ), "ProtoMS tools file (%s) is not present in the expected place." %tools_files)

#Entry point for nosetests or unittests


if  __name__ == '__main__':
    unittest.main()
    nose.runmodule()
