import unittest

import framework


class Jaws1Test(framework.BaseTest):
    ref_dir = "tests/jaws1/"

    copy_files = [
        "fragment.pdb",
        "protein.pdb",
        "water.pdb"
    ]

    setup_args = [
        "-s", "jaws1",
        "-l", "fragment.pdb",
        "-p", "protein.pdb",
        "--nequil", "0",
        "--nprod", "100",
        "--ranseed", "100001",
        "--setupseed", "100000",
        "--dumpfreq", "10",
        "-w water.pdb"
    ]

    setup_output_files = [
        "water_clr.pdb",
        "jaws1_box.pdb",
        "jaws1_wat.pdb",
        "run_jaws.cmd"
    ]

    simulation_args = [
        "run_jaws.cmd"
    ]

    simulation_output_directory = "out_jaws"

    simulation_output_files = [
        "results",
        "accept",
        "all.pdb",
        "restart",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
