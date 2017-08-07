"""
Short script to remove the TER line separating protein chains
This allows ProtoMS to read multi-chain proteins as a single long chain
The discontinuities between chain ends/beginnings are not an issue as ProtoMS will assume this is part of a scoop

Marley Samways
"""


def remove_ter(pdb_in, pdb_out):
    '''
    Remove the TER lines separating protein chains
    '''
    # Read input
    in_file = open(pdb_in, 'r')
    pdb_lines = in_file.readlines()
    in_file.close()
    # Scan file to remove TER lines (only those separating protein chains)
    ter = True
    while ter:  # Keep searching whilst TER is true...
        found = False  # Marks whether a suitable TER line has been found
        for i in range(1, len(pdb_lines)-1):
            # Check that this is a TER line
            curr_line = pdb_lines[i]
            if curr_line[:3] != "TER":
                continue
            # Need to check if previous and following lines are protein lines...
            prev_line = pdb_lines[i-1]
            next_line = pdb_lines[i+1]
            prev_resname = prev_line[17:20]
            next_resname = next_line[17:20]
            protein_resnames = ['GLH', 'ILE', 'GLN', 'GLY', 'GLU', 'HIP', 'HIS', 'SER', 'LYS', 'PRO', 'CYX', 'HIE', 'LYN', 'ASH', 'ASN', 'CYS', 'VAL', 'THR', 'ASP', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
            if prev_resname in protein_resnames and next_resname in protein_resnames:
                 found = True
                 ter_i = i
                 break
            else:
                 continue
        if found:
            pdb_lines = pdb_lines[:i] + pdb_lines[i+1:]
        else:
            ter = False
    # After removing all TER lines separating protein chains, the data must be written out
    out_file = open(pdb_out, 'w')
    for line in pdb_lines:
        out_file.write(line)
    out_file.close()


def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="Script to ensure that no TER lines are separating protein lines")
    parser.add_argument('-i', '--input', help="Input PDB file")
    parser.add_argument('-o', '--output', help="Output PDB file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    remove_ter(args.input, args.output)




