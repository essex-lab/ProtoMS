import unittest

import framework


class Jaws1SetupTest(framework.BaseTest):
    ref_dir = "tests/jaws1/"

    input_files = [
        "fragment.pdb",
        "protein.pdb",
        "water.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "jaws1",
        "-l", "fragment.pdb",
        "-p", "protein.pdb",
        "--nequil", "0",
        "--nprod", "100",
        "--ranseed", "100001",
        "--setupseed", "100000",
        "--dumpfreq", "10",
        "-w", "water.pdb",
        "--gaff", "gaff14"
    ]

    output_files = [
        "water_clr.pdb",
        "jaws1_box.pdb",
        "jaws1_wat.pdb",
        "run_jaws.cmd",
        "fragment_box.pdb",
        "fragment.frcmod",
        "fragment.prepi",
        "fragment.tem",
        "fragment.zmat",
        "protein_scoop.pdb"
    ]


class Jaws1SimulationTest(framework.BaseTest):
    ref_dir = "tests/jaws1/"

    input_files = [
        "fragment.pdb",
        "water_clr.pdb",
        "jaws1_box.pdb",
        "jaws1_wat.pdb",
        "run_jaws.cmd",
        "fragment.tem",
        "protein_scoop.pdb"
    ]

    executable = "build/protoms3"

    args = [
        "run_jaws.cmd"
    ]

    output_directories = [
        "out_jaws"
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
