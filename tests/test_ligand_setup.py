# This test is to check the protoms.py ligand setup.

import nose
import unittest
import os
import site

# ---------------------------------------------
# Ligand Setup test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]

site.addsitedir(proto_env)
from protoms import _get_prefix, _load_ligand_pdb, _prep_ligand
from tools import simulationobjects


class MockArgs:
    template = None


class TestLigSetup(unittest.TestCase):
    def setUp(self):
        super(TestLigSetup, self).setUp()
        self.tarlist = []   # List of files to be stored
        self.prefix = None  # Ligand filename (prefix)
        self.ligands = []   # List of ligands
        self.ligpdb = None  # List of ligand pdb-files
        self.ligtem = None  # List of ligand template files
        self.ligobj = None  # Merged pdb object of all ligand pdb objects
        self.ligfiles = {}  # Dictionary of filenames for ligand

    def tearDown(self):
        super(TestLigSetup, self).tearDown()

    def test_prep_ligand(self):
        try:
            self.prefix = _get_prefix("dcb.pdb")
            self.ligands.append(self.prefix)

            self.ligfiles[self.prefix] = {}

            self.ligfiles[self.prefix]["pdb"], self.ligfiles[self.prefix]["obj"] = _load_ligand_pdb(self.prefix, [os.path.join(proto_env, "tests/setup")])

            self.ligobj = self.ligfiles[self.ligands[0]]["obj"]  # Unmerged pdb object for single ligand

            args = MockArgs()
            _prep_ligand(self.ligfiles[self.prefix], 0, 0, self.ligobj, [" "], self.tarlist, args)

            print(self.tarlist)

        except ImportError as e:
            print(e)
            print("Ligand file not found.")


if __name__ == '__main__':
    logger = simulationobjects.setup_logger('protoms_py.log')
    unittest.main()
    nose.runmodule()
