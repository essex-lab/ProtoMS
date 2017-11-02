"""This script contains tests developed during implementation of an
atom-by-atom softcore functionality. Unilke other tests for ProtoMS,
the below do not compare outputs against set reference data but look
for internal consistency within a set of simulations.

All tests are performed with a simple idealised system as shown below,
containing a t3p water and a bromine dimer aligned within a plane.

              H
               \                Br
                \               |
                 O              |
                /               |
               /                Br
              H

Three energy calculations are performed at lambda=0.5, with standard
softcore parameters on LJ and Coulomb interactions:
a) Both Br atoms have softcores applied
b) Neither Br atoms have softcores applied
c) One atom has softcores applied

The formulation of the test means that, if correctly implemented, the
interaction energies of the three calculations should obey the below
rule.

(a+b)/2 = c

Letters here correspond to entries in the list of calculations
above. Testing this relationship is the basis of this
script. Individual tests for the interaction of solutes with each
molecule type in ProtoMS are required.

"""

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


class EnergiesSimulationSoftProteinTest(EnergiesSimulationSoftSolventTest):
    """In order to test the protein-solute energy code, the t3p water
    model is implemented as a protein residue in t3p.ff. This required
    the inclusion of additional 'backbone' atoms that are assigned
    dummy parameters so do not influence the calculation of energies.
    """

    input_files = [
        "run_soft_protein.cmd",
        "run_mixed_protein.cmd",
        "run_protein.cmd",
        "soft_solute.pdb",
        "soft_protein.pdb",
        "brd.tem",
        "t3p.ff"
    ]

    args = [
        ["run_soft_protein.cmd"],
        ["run_mixed_protein.cmd"],
        ["run_protein.cmd"]
    ]

    output_directories = [
        "out_soft_protein",
        "out_mixed_protein",
        "out_protein"
    ]

    def get_result_energies(self, filename):
        snapshot = sim.SnapshotResults()
        with open(filename) as f:
            snapshot.parse(f)
        term = snapshot.interaction_energies['protein1-brd1']
        return (term[0].curr, term[1].curr)


if __name__ == '__main__':
    unittest.main()
