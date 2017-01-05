import unittest

import framework


class ProtSetupTest(framework.BaseTest):
    ref_dir = "tests/setup/"

    copy_files = [
        "dcb.pdb",
        "protein.pdb"
    ]

    setup_args = [
        "-s", "none",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--setupseed", "100000",
        "--gaff", "gaff14"
    ]

    setup_output_files = [
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
