import unittest

import framework


class EnergiesSimulationTip3pTest(framework.BaseTest):
    ref_dir = "tests/energies/"

    input_files = [
        "run_t3p.cmd",
        "protein_scoop.pdb",
        "t3p.pdb"
    ]

    executable = "build/protoms3"

    args = [
        "run_t3p.cmd"
    ]

    output_directories = [
        "out_t3p"
    ]

    output_files = [
        "info",
        "results",
        "warning"
    ]


class EnergiesSimulationTip4pTest(framework.BaseTest):
    ref_dir = "tests/energies/"

    input_files = [
        "run_t4p.cmd",
        "protein_scoop.pdb",
        "t4p.pdb"
    ]

    executable = "build/protoms3"

    args = [
        "run_t4p.cmd"
    ]

    output_directories = [
        "out_t4p"
    ]

    output_files = [
        "info",
        "results",
        "warning"
    ]


if __name__ == '__main__':
    unittest.main()
