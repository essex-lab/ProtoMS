import os
import sys
import unittest

import framework

sys.path.append(os.path.join(os.environ['PROTOMSHOME'], 'tools'))
import simulationobjects as sim


class EnergiesSimulationSoftSolventTest(framework.BaseTest):
    ref_dir = "tests/energies/"

    input_files = [
        "run_soft_solvent.cmd",
        "run_mixed_solvent.cmd",
        "run_solvent.cmd",
        "soft_solute.pdb",
        "soft_solvent.pdb",
        "brd.tem"
    ]

    executable = "build/protoms3"

    args = [
        ["run_soft_solvent.cmd"],
        ["run_mixed_solvent.cmd"],
        ["run_solvent.cmd"]
    ]

    output_directories = [
        "out_soft_solvent",
        "out_mixed_solvent",
        "out_solvent"
    ]

    output_files = [
        "results",
    ]

    def test(self):
        self._helper_copy_input_files()

        print("\nTEST_RUN\n")
        for args in self.args:
            args = [self.executable] + args
            if self.mpi_processes > 0:
                args = ["mpirun", "-np", str(self.mpi_processes)] + args
            self._helper_subprocess_call(args)

        print("\nTEST_OUTPUT\n")
        out_files = [os.path.join(d, f)
                     for f in self.output_files
                     for d in self.output_directories]
        soft_energies = self.get_result_energies(out_files[0])
        mixed_energies = self.get_result_energies(out_files[1])
        energies = self.get_result_energies(out_files[2])
        for se, me, e in zip(soft_energies, mixed_energies, energies):
            self.assertAlmostEqual((se+e)/2, me, 2)

    def get_result_energies(self, filename):
        snapshot = sim.SnapshotResults()
        with open(filename) as f:
            snapshot.parse(f)
        term = snapshot.interaction_energies['brd-solvent']
        return (term[0].curr, term[1].curr)


if __name__ == '__main__':
    unittest.main()
