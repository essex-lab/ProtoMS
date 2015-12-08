"""This test is to check the ProtoMS setup of simulations with multiple lambda windows and the MPI implementation of protoms3 program."""

import nose
import unittest
import argparse
import os
import sys
import subprocess
import logging
import time
import re
import numpy as np
import filecmp
import protoms

from protoms import _is_float, _get_prefix, _locate_file, _merge_templates
from protoms import _load_ligand_pdb, _prep_ligand, _prep_protein, _prep_singletopology
from protoms import _prep_gcmc, _prep_jaws2, _cleanup, _wizard

import tools
from tools import simulationobjects

from subprocess import call

# --------------------------------------------------------
# ProtoMS setup and RETI single topology simulations test
# --------------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/RETI_sngl/"

output_files_setup = ["ethane_box.pdb", "ethtmeo_ele.tem", "ethtmeo_vdw.tem", "ethtmeo_comb.tem"]
outfiles_setup = ["run_comb_free.cmd", "run_comb_gas.cmd", "run_ele_free.cmd", "run_ele_gas.cmd", "run_vdw_free.cmd", "run_vdw_gas.cmd"]
out_sim_files = ["results", "accept", "all.pdb", "restart.prev", "warning", "info", "restart", "results_inst"]


class TestRETIsngl(unittest.TestCase):
    """Test for RETI/MPI function."""

    def setUp(self):
        super(TestRETIsngl, self).setUp()

    def tearDown(self):
        super(TestRETIsngl, self).tearDown()

    def test_RETI_sngl(self):
        """Test for RETI single topology/MPI function."""

        if((call("python2.7 $PROTOMSHOME/protoms.py -s singletopology -l ethane.pdb methanol.pdb --nequil 0 --nprod 10 --lambdas 0.00 0.33 0.67 1.00 --ranseed 100000 --dumpfreq 1 --cleanup --singlemap single_cmap.dat", shell=True)) == 0):
            # Checking whether the required output files have been setup for RETI single topology protoms.py setup.

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            for outfile in outfiles_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            # Checking content of RETI single topology setup output files with reference data files.
            for outfile in output_files_setup:
                self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                "Content mismatch between output and reference for file {0}".format(outfile))

            for out_files in outfiles_setup:
                if((call("bash content_cmd_comp.sh", shell=True)) == 0):
                    print("\n Content matched for %s." % out_files)
                else:
                    raise ValueError("Content mismatch between output and reference %s" % (out_files))

        else:
            raise simulationobjects.SetupError("ProtoMS setup for RETI single topology and command files generation failed.")

        if((call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_comb_free.cmd", shell=True)) == 0):
            # Checking whether the simulation output files have been created successfully for RETI free phase leg of a single topology for combined perturbation.

            if(os.path.exists("out_comb_free")):
                for root, dirs, files in os.walk("out_comb_free"):
                    if len(dirs) != 0:
                        for d in dirs:
                            for out_files in out_sim_files:
                                outfile_rel = os.path.join("out_comb_free", d, out_files)
                                self.assertTrue(os.path.exists(outfile_rel),
                                                "Simulation file {0} is missing.".format(outfile_rel))

                                # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                                if out_files == "info":
                                    if((call("bash content_info_freecomp.sh", shell=True)) == 0):
                                        print("\n Free phase leg info files content matched.")
                                    else:
                                        raise ValueError("Content mismatch between output and reference info file for lambda value {0}.".format(d))
                                else:
                                    self.assertTrue(filecmp.cmp(outfile_rel, os.path.join(ref_dir, outfile_rel)),
                                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("RETI free phase leg of a single topology simulation is not successful.")

        if((call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_comb_gas.cmd", shell=True)) == 0):
            # Checking whether the simulation output files have been created successfully for RETI gas phase leg of a single topology for combined perturbation.
            if(os.path.exists("out_comb_gas")):
                for root, dirs, files in os.walk("out_comb_gas"):
                    if len(dirs) != 0:
                        for d in dirs:
                            for out_files in out_sim_files:
                                outfile_rel = os.path.join("out_comb_gas", d, out_files)
                                self.assertTrue(os.path.exists(outfile_rel),
                                                "Simulation file {0} is missing.".format(outfile_rel))

                                # Checking content of RETI gas phase leg of a single topology simulation output files with reference data.
                                if out_files == "info":
                                    if((call("bash content_info_gascomp.sh", shell=True)) == 0):
                                        print("\n Gas phase leg info files matched.")
                                        continue
                                    else:
                                        raise ValueError("Content mismatch between output and reference info file for lambda value {0}.".format(d))
                                else:
                                    self.assertTrue(filecmp.cmp(outfile_rel, os.path.join(ref_dir, outfile_rel)),
                                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("RETI gas phase leg of a single topology simulation was not successful.")


# Entry point to nosetests or unittests
if __name__ == "__main__":
    unittest.main()
    nose.runmodule()
