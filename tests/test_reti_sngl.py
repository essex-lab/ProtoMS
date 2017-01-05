import unittest

import framework


class RetiSnglSetupTest(framework.BaseTest):
    ref_dir = "tests/RETI_sngl/"

    copy_files = [
        "ethane.pdb",
        "methanol.pdb",
        "single_cmap.dat"
    ]

    setup_args = [
        "-s", "singletopology",
        "-l", "ethane.pdb", "methanol.pdb",
        "--lambdas", "0.00", "0.33", "0.67", "1.00",
        "--nequil", "0",
        "--nprod", "10",
        "--ranseed", "100000",
        "--dumpfreq", "1",
        "--cleanup",
        "--singlemap", "single_cmap.dat",
        "--gaff", "gaff14"
    ]

    setup_output_files = [
        "ethane_box.pdb",
        "ethtmeo_comb.tem",
        "ethtmeo_ele.tem",
        "ethtmeo_vdw.tem",
        "run_comb_free.cmd",
        "run_comb_gas.cmd",
        "run_ele_free.cmd",
        "run_ele_gas.cmd",
        "run_vdw_free.cmd",
        "run_vdw_gas.cmd"
    ]


class RetiSnglSimulationFreeTest(framework.BaseTest):
    ref_dir = "tests/RETI_sngl/"

    copy_files = [
        "run_comb_free.cmd",
        "ethane.pdb",
        "ethtmeo_comb.tem",
        "ethane_box.pdb"
    ]

    simulation_mpi_processes = 4

    simulation_args = [
        "run_comb_free.cmd"
    ]

    simulation_output_directories = [
        "out_comb_free/lam-0.000",
        "out_comb_free/lam-0.330",
        "out_comb_free/lam-0.670",
        "out_comb_free/lam-1.000"
    ]

    simulation_output_files = [
        "accept",
        "all.pdb",
        "info",
        "restart",
        "restart.prev",
        "results",
        "results_inst",
        "warning"
    ]


class RetiSnglSimulationGasTest(framework.BaseTest):
    ref_dir = "tests/RETI_sngl/"

    copy_files = [
        "run_comb_gas.cmd",
        "ethane.pdb",
        "ethtmeo_comb.tem"
    ]

    simulation_mpi_processes = 4

    simulation_args = [
        "run_comb_gas.cmd"
    ]

    simulation_output_directories = [
        "out_comb_gas/lam-0.000",
        "out_comb_gas/lam-0.330",
        "out_comb_gas/lam-0.670",
        "out_comb_gas/lam-1.000"
    ]

    simulation_output_files = [
        "accept",
        "all.pdb",
        "info",
        "restart",
        "restart.prev",
        "results",
        "results_inst",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
