import unittest

import framework


class RetiDblSetupTest(framework.BaseTest):
    ref_dir = "tests/RETI_dbl/"

    input_files = [
        "ethane.pdb",
        "methanol.pdb"
    ]

    executable = "protoms.py"

    args = [
        "-s", "dualtopology",
        "-l", "ethane.pdb", "methanol.pdb",
        "--lambdas", "0.00", "0.33", "0.67", "1.00",
        "--nequil", "0",
        "--nprod", "10",
        "--ranseed", "100000",
        "--dumpfreq", "1",
        "--cleanup",
        "--gaff", "gaff14",
    ]

    output_files = [
        "ethane_box.pdb",
        "eth-meo.tem",
        "run_free.cmd"
    ]


class RetiDblSimulationTest(framework.BaseTest):
    ref_dir = "tests/RETI_dbl/"

    input_files = [
        "run_free.cmd",
        "eth-meo.tem",
        "ethane.pdb",
        "methanol.pdb",
        "ethane_box.pdb"
    ]

    mpi_processes = 4

    executable = "build/protoms3"

    args = [
        "run_free.cmd"
    ]

    output_directories = [
        "out_free/lam-0.000",
        "out_free/lam-0.330",
        "out_free/lam-0.670",
        "out_free/lam-1.000"
    ]

    output_files = [
        "accept",
        "all.pdb",
        "restart",
        "restart.prev",
        "results",
        "results_inst",
        "warning",
        "info"
    ]


if __name__ == '__main__':
    unittest.main()
