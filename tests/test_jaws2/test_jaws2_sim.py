"""This test is to check the ProtoMS setup and protoms3 run of JAWS stage 2 simulations."""

import nose
import unittest
import os
import filecmp
import site

from subprocess import call

# ------------------------------------------------
# ProtoMS setup and JAWS stage 2 simulations test
# ------------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/jaws2/"

site.addsitedir(proto_env)
from tools import simulationobjects
from tools import testtools

output_files_setup = ["fragment.tem", "fragment.frcmod", "fragment.prepi", "fragment.zmat", "fragment_box.pdb", "protein_scoop.pdb"]
outfiles_setup = ["jaws2_wat1.pdb", "jaws2_wat2.pdb", "jaws2_wat3.pdb", "jaws2_not1.pdb", "jaws2_not2.pdb", "jaws2_not3.pdb", "run_jaws2-w1_jaws.cmd", "run_jaws2-w2_jaws.cmd", "run_jaws2-w3_jaws.cmd"]
output_subdirs = ["j_+8.000", "j_+10.000", "j_+12.000", "j_+14.000"]
out_sim_files = ["accept", "all.pdb", "info", "restart", "restart.prev", "results", "warning"]


class TestJAWS2(unittest.TestCase):
    """Test for JAWS2 function."""

    def setUp(self):
        super(TestJAWS2, self).setUp()

    def tearDown(self):
        super(TestJAWS2, self).tearDown()

    def test_jaws2(self):
        """Test for JAWS2 function."""
        comparetools = testtools.CompareTools(ref_dir, verbose=True)

        if((call("python2.7 $PROTOMSHOME/protoms.py -s jaws2 -l fragment.pdb -p protein.pdb --gcmcwater jaws2_waters.pdb --jawsbias 8 10 12 14 --nequil 0 --nprod 100 --ranseed 100000 --setupseed 100000 --dumpfreq 10", shell=True)) == 0):
            # Checking whether the required output files have been setup for JAWS Stage 2 protoms.py setup.

            for outfile in output_files_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            for outfile in outfiles_setup:
                self.assertTrue(os.path.exists(outfile),
                                "ProtoMS setup output file {0} is missing.".format(outfile))

            # Checking content of setup output and reference files for JAWS Stage II.
            for outfile in output_files_setup:
                if outfile == "protein_scoop.pdb":
                    if call("bash content_ps_comp.sh", shell=True) == 0:
                        print("\n Protein scoop files matched.")
                    else:
                        raise ValueError("Content mismatch between output and reference Protein scoop files.")
                else:
                    self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile))

            for outfile in outfiles_setup:
                if outfile == "run_bnd.cmd":
                    if call("bash content_cmd_comp.sh", shell=True) == 0:
                        print("\n Output Command files matched.")
                    else:
                        raise ValueError("Content mismatch between output and reference command files.")
                else:
                    self.assertTrue(filecmp.cmp(outfile, os.path.join(ref_dir, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile))

        else:
            raise simulationobjects.SetupError("ProtoMS setup and command files generation for JAWS Stage 2 failed.")

        if call("mpirun -np 4 $PROTOMSHOME/build/protoms3 run_jaws2-w1_jaws.cmd", shell=True) == 0:
            # Checking whether the simulation output files have been created successfully for JAWS Stage 2.
            for subdir in output_subdirs:
                for outfile in out_sim_files:
                    outfile_rel = os.path.join("out_jaws2-w1", subdir, outfile)
                    self.assertTrue(os.path.exists(outfile_rel),
                                    "Simulation file {0} is missing.".format(outfile_rel))

                    # Checking content of RETI free phase leg of a single topology simulation output files with reference data.
                    self.assertTrue(comparetools.compare((outfile_rel, outfile)),
                                    "Content mismatch between output and reference for file {0}".format(outfile_rel))

        else:
            raise simulationobjects.SetupError("JAWS Stage 2 check simulation was not successful.")

# Entry point to unittests or nosetests
if __name__ == "__main__":
    unittest.main()
    nose.runmodule()
