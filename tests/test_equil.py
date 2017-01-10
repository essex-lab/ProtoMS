import unittest

import framework


class EquilSetupTest(framework.BaseTest):
    ref_dir = "tests/equil/"

    input_files = [
        "dcb.frcmod",
        "dcb.pdb",
        "dcb.prepi",
        "dcb.tem",
        "dcb.zmat",
        "dcb_box.pdb",
        "protein.pdb",
        "protein_scoop.pdb",
        "water.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "equilibration",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--nequil", "100",
        "--ranseed", "100000",
        "--gaff", "gaff14"
    ]

    output_files = [
        "run_bnd.cmd"
    ]


class EquilSimulationTest(framework.BaseTest):
    ref_dir = "tests/equil/"

    input_files = [
        "dcb.frcmod",
        "dcb.pdb",
        "dcb.prepi",
        "dcb.tem",
        "dcb.zmat",
        "dcb_box.pdb",
        "protein.pdb",
        "protein_scoop.pdb",
        "water.pdb",
        "run_bnd.cmd"
    ]

    executable = "build/protoms3"

    args = [
        "run_bnd.cmd"
    ]

    output_directories = [
        "out_bnd"
    ]

    output_files = [
        "equil_bnd.pdb",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
