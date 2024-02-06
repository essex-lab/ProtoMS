from protomslib import simulationobjects as sim
import argparse

parser = argparse.ArgumentParser(
    description="Script to modify ProtoMS simulation output "
    "trajectories produced using GCMC such that a constant "
    "number of atoms is present in every frame. The resulting pdb "
    "can then be viewed with VMD."
)
parser.add_argument(
    "-i", "--in_file", required=True, help="pdb trajectory input file"
)
parser.add_argument(
    "-o", "--out_file", required=True, help="pdb trajectory output file"
)
parser.add_argument(
    "-n",
    "--nwat",
    required=True,
    help="Largest number of inserted " " GC waters in input trajectory file",
    type=int,
)
args = parser.parse_args()


def make_ref(index):
    """Create a new reference water whenever one is
    missing from the GCMC trajectory."""
    ref_res = sim.Residue(name="WA1", index=index)
    ref_res.addAtom(
        sim.Atom(
            index=1,
            name="O00",
            resname="WA1",
            resindex=index,
            coords=[1000.0 + index * 3, 0.0, 0.0],
        )
    )
    ref_res.addAtom(
        sim.Atom(
            index=1,
            name="H01",
            resname="WA1",
            resindex=index,
            coords=[1000.0 + index * 3, 1.0, 0.0],
        )
    )
    ref_res.addAtom(
        sim.Atom(
            index=1,
            name="H02",
            resname="WA1",
            resindex=index,
            coords=[1000.0 + index * 3, 0.0, 1.0],
        )
    )
    ref_res.addAtom(
        sim.Atom(
            index=1,
            name="M03",
            resname="WA1",
            resindex=index,
            coords=[1000.0 + index * 3, 0.1, 0.1],
        )
    )
    return ref_res


s = sim.PDBSet()
s.read(args.in_file)

for pdb in s.pdbs:
    count = len([r for r in pdb.residues if pdb.residues[r].name == "WA1"])

    for i in range(args.nwat - count):
        try:
            index = sorted(pdb.residues.keys())[-1] + 1
        except IndexError:
            index = 1
        pdb.residues[index] = make_ref(index)

s.write(args.out_file)
