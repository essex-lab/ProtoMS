import unittest

import framework


class GcmcSetupBoxTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    copy_files = [
        "wat.pdb"
    ]

    setup_executable = "/home/james/projects/protoms-dev/tools/make_gcmcbox.py"

    setup_args = [
        "-s", "wat.pdb",
    ]

    setup_output_files = [
        "gcmc_box.pdb",
    ]


class GcmcSetupTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    copy_files = [
        "protein.pdb",
        "water.pdb",
        "wat.pdb",
        "gcmc_box.pdb"
    ]

    setup_args = [
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

    setup_output_files = [
        "gcmc_wat.pdb",
        "run_gcmc.cmd"
    ]


class GcmcSimulationTest(framework.BaseTest):
    ref_dir = "tests/gcmc/"

    copy_files = [
        "run_gcmc.cmd",
        "protein.pdb",
        "gcmc_wat.pdb",
        "water.pdb"
    ]

    simulation_mpi_processes = 2

    simulation_args = [
        "run_gcmc.cmd"
    ]

    simulation_output_directories = [
        "out_gcmc/b_+19.000",
        "out_gcmc/b_+20.000"
    ]

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
