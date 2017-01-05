import unittest

import framework


class EquilSetupTest(framework.BaseTest):
    ref_dir = "tests/equil/"

    copy_files = [
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

    setup_args = [
        "-s", "equilibration",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--nequil", "100",
        "--ranseed", "100000",
        "--gaff", "gaff14"
    ]

    setup_output_files = [
        "run_bnd.cmd"
    ]


class EquilSimulationTest(framework.BaseTest):
    ref_dir = "tests/equil/"

    copy_files = [
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

    simulation_args = [
        "run_bnd.cmd"
    ]

    simulation_output_directory = "out_bnd"

    simulation_output_files = [
        "equil_bnd.pdb",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
