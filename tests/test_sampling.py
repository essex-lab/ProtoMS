import unittest

import framework


class SamplingSetupTest(framework.BaseTest):
    ref_dir = "tests/sampling/"

    copy_files = [
        "dcb.pdb",
        "protein.pdb",
        "water.pdb"
    ]

    setup_args = [
        "-s", "sampling",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--nequil", "0",
        "--nprod", "100",
        "--ranseed", "100000",
        "--dumpfreq", "10",
        "--gaff", "gaff14"
    ]

    setup_output_files = [
        "dcb.frcmod",
        "dcb.prepi",
        "dcb.tem",
        "dcb.zmat",
        "dcb_box.pdb",
        "protein_scoop.pdb",
        "run_bnd.cmd"
    ]


class SamplingSimulationTest(framework.BaseTest):
    ref_dir = "tests/sampling/"

    copy_files = [
        "dcb.pdb",
        "dcb.tem",
        "protein_scoop.pdb",
        "water.pdb",
        "run_bnd.cmd"
    ]

    simulation_args = [
        "run_bnd.cmd"
    ]

    simulation_output_directory = "out_bnd"

    simulation_output_files = [
        "accept",
        "all.pdb",
        "info",
        "restart",
        "restart.prev",
        "results",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
