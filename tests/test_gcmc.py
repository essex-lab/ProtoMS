import unittest

import framework


class GcmcSetupBoxTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    input_files = [
        "wat.pdb"
    ]

    executable = "tools/make_gcmcbox.py"

    args = [
        "-s", "wat.pdb",
    ]

    output_files = [
        "gcmc_box.pdb",
    ]


class GcmcSetupTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    input_files = [
        "protein.pdb",
        "water.pdb",
        "wat.pdb",
        "gcmc_box.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "gcmc",
        "-sc", "protein.pdb",
        "-w", "water.pdb",
        "--gcmcwater", "wat.pdb",
        "--gcmcbox", "gcmc_box.pdb",
        "--adams", "19", "20",
        "--nequil", "0",
        "--nprod", "1000",
        "--ranseed", "100000",
        "--dumpfreq", "100",
        "--capradius", "26",
        "--gaff", "gaff14"
    ]

    output_files = [
        "gcmc_wat.pdb",
        "run_gcmc.cmd"
    ]


class GcmcSimulationTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    input_files = [
        "run_gcmc.cmd",
        "protein.pdb",
        "gcmc_wat.pdb",
        "water.pdb"
    ]

    mpi_processes = 2

    executable = "build/protoms3"

    args = [
        "run_gcmc.cmd"
    ]

    output_directories = [
        "out_gcmc/b_+19.000",
        "out_gcmc/b_+20.000"
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
