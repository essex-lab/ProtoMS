import unittest

import framework


class EnergiesSimulationTip3pTest(framework.BaseTest):
    ref_dir = "tests/energies/"

    copy_files = [
        "run_t3p.cmd",
        "protein_scoop.pdb",
        "t3p.pdb"
    ]

    simulation_args = [
        "run_t3p.cmd"
    ]

    simulation_output_directory = "out_t3p"

    simulation_output_files = [
        "info",
        "results",
        "warning"
    ]


class EnergiesSimulationTip4pTest(framework.BaseTest):
    ref_dir = "tests/energies/"

    copy_files = [
        "run_t4p.cmd",
        "protein_scoop.pdb",
        "t4p.pdb"
    ]

    simulation_args = [
        "run_t4p.cmd"
    ]

    simulation_output_directory = "out_t4p"

    simulation_output_files = [
        "info",
        "results",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
