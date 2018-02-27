import unittest

import framework


class ProtSetupTest(framework.BaseTest):
    ref_dir = "tests/setup/"

    input_files = [
        "dcb.pdb",
        "protein.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "none",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--charge", "0",
        "--setupseed", "100000",
        "--gaff", "gaff14"
    ]

    output_files = [
        "dcb.prepi",
        "dcb.frcmod",
        "dcb.zmat",
        "dcb.tem",
        "dcb_box.pdb",
        "protein_scoop.pdb",
        "water.pdb"
    ]


if __name__ == '__main__':
    unittest.main()
