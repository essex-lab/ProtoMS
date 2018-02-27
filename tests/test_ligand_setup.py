# This test is to check the protoms.py ligand setup.

import unittest

import framework


class LigandSetupTest(framework.BaseTest):
    ref_dir = "tests/setup"

    input_files = ['dcb.pdb']

    executable = "protoms.py"

    args = ['-l', 'dcb.pdb',
            '--charge', '0',
            '--gaff', 'gaff14',
            '--setupseed', "100000"]

    output_files = ['dcb.prepi', 'dcb.frcmod', 'dcb.zmat',
                    'dcb.tem', 'dcb_box.pdb']


if __name__ == '__main__':
    unittest.main()
