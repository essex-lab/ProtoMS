#!/usr/bin/env python

"""
Tests for ProtoMS Tools (protoms.py): Setup tests
"""

import os
import subprocess
import unittest
import nose

# ---------------------------------------------
# ProtoMS tools tests
# ---------------------------------------------

# List containing ProtoMS tools.
protoms_tools = [
    "ambertools.py",
    "build_template.py",
    "calc_clusters.py",
    "calc_density.py",
    "calc_dg_cycle.py",
    "calc_dg.py",
    "calc_gcap_surface.py",
    "calc_gci_reweight.py",
    "calc_gci.py",
    "calc_replicapath.py",
    "calc_rmsd.py",
    "calc_series.py",
    "calc_ti_decomposed.py",
    "clear_gcmcbox.py",
    "convertatomnames.py",
    "convertwater.py",
    "distribute_waters.py",
    "divide_pdb.py",
    "generate_input.py",
    "make_dummy.py",
    "make_gcmcbox.py",
    "make_gcmc_traj.py",
    "make_single.py",
    "merge_templates.py",
    "plot_theta.py",
    "scoop.py",
    "solvate.py",
    "split_jawswater.py",
]


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
            path = os.environ["PROTOMSHOME"] + "/tools/" + tools_files
            self.assertTrue(
                os.path.exists(path),
                "ProtoMS tools file {0} is not present in the expected place.".format(
                    tools_files
                ),
            )
            # at the very least test that we can call each script
            # and get a help string
            with open(os.devnull, "w") as f:
                subprocess.check_call(["python", path, "--help"], stdout=f)


# Entry point for nosetests or unittests
if __name__ == "__main__":
    unittest.main()
    nose.runmodule()
