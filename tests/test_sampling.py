import unittest

import framework


class SamplingSetupTest(framework.BaseTest):
    ref_dir = "tests/sampling/"

    input_files = [
        "dcb.pdb",
        "protein.pdb",
        "water.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "sampling",
        "-l", "dcb.pdb",
        "-p", "protein.pdb",
        "--nequil", "0",
        "--nprod", "100",
        "--ranseed", "100000",
        "--dumpfreq", "10",
        "--gaff", "gaff14"
    ]

    output_files = [
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

    input_files = [
        "dcb.pdb",
        "dcb.tem",
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
