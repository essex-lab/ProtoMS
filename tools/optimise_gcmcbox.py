"""
Optimise rotations of a ligand to minimise volume of GCMC box
File needs some tidying...

Marley Samways
"""


import numpy as np
import logging


def optimise_system(ligand_file, other_files=None, heavy=False):
    """
    Sample rotational angles to identify those which give the smallest box around the ligand
    Ligand and protein will then be rotated by these angles and PDB files will be written
    Currently not written in the simulationobjects framework - may change this later...
    """
    min_volume = box_volume(ligand_file, heavy=heavy)
    opt_angles = np.zeros(3)
    for x in range(0, 370, 10):
        for y in range(0, 370, 10):
            for z in range(0, 190, 10):
                angles = np.array([float(x), float(y), float(z)])
                angles *= np.pi / 180.0
                # print angles
                test_volume = box_volume(ligand_file, angles=angles, heavy=heavy)
                # print test_volume
                if test_volume < min_volume:
                    min_volume = test_volume
                    # print opt_angles
                    opt_angles = angles
    # Need to optimise all files and rename original input
    input_files = [ligand_file]
    if other_files != None:
        for file in other_files:
            input_files.append(file)
    orig_files = []
    for file in input_files:
        # Need to rename the original PDB files, such that they are not lost..
        new_file = file[:-4] + '_orig.pdb'
        orig_files.append(new_file)
    for i in range(len(input_files)):
        # Making a copy of the input
        write_pdb(input_files[i], orig_files[i], angles=np.zeros(3))
        # Overwriting the original input with the rotated structure
        write_pdb(input_files[i], input_files[i], angles=opt_angles)


def rotate_coords(atom_coords, angles, forward=True):
    """
    Rotate a set of coordinates by a set of angles.
    If forward is set to false, the rotation is set to work in reverse
    """
    if forward:
        new_coords = np.zeros(3)
        for i in range(3):
            new_coords[i] = atom_coords[i]
        # Rotate about x
        new_coords[1] = (atom_coords[1]*np.cos(angles[0])) - (atom_coords[2]*np.sin(angles[0]))
        new_coords[2] = (atom_coords[1]*np.sin(angles[0])) + (atom_coords[2]*np.cos(angles[0]))
        for i in range(3):
            atom_coords[i] = new_coords[i]
        # Rotate about y
        new_coords[0] = (atom_coords[0]*np.cos(angles[1])) + (atom_coords[2]*np.sin(angles[1]))
        new_coords[1] = atom_coords[1]
        new_coords[2] = (atom_coords[2]*np.cos(angles[1])) - (atom_coords[0]*np.sin(angles[1]))
        for i in range(3):
            atom_coords[i] = new_coords[i]
        # Rotate about z
        new_coords[0] = (atom_coords[0]*np.cos(angles[2])) - (atom_coords[1]*np.sin(angles[2]))
        new_coords[1] = (atom_coords[0]*np.sin(angles[2])) + (atom_coords[1]*np.cos(angles[2]))
        new_coords[2] = atom_coords[2]
    else:
        new_coords = np.zeros(3)
        for i in range(3):
            new_coords[i] = atom_coords[i]
        # Rotate about z
        new_coords[0] = (atom_coords[0]*np.cos(-angles[2])) - (atom_coords[1]*np.sin(-angles[2]))
        new_coords[1] = (atom_coords[0]*np.sin(-angles[2])) + (atom_coords[1]*np.cos(-angles[2]))
        new_coords[2] = atom_coords[2]
        for i in range(3):
            atom_coords[i] = new_coords[i]
        # Rotate about y
        new_coords[0] = (atom_coords[0]*np.cos(-angles[1])) + (atom_coords[2]*np.sin(-angles[1]))
        new_coords[1] = atom_coords[1]
        new_coords[2] = (atom_coords[2]*np.cos(-angles[1])) - (atom_coords[0]*np.sin(-angles[1]))
        for i in range(3):
            atom_coords[i] = new_coords[i]
        # Rotate about x
        new_coords[0] = atom_coords[0]
        new_coords[1] = (atom_coords[1]*np.cos(-angles[0])) - (atom_coords[2]*np.sin(-angles[0]))
        new_coords[2] = (atom_coords[1]*np.sin(-angles[0])) + (atom_coords[2]*np.cos(-angles[0]))
    return new_coords


def box_volume(pdb_file, angles=None, heavy=False):
    """
    Figure out volume of around the ligand, given some rotation
    """
    min_coords = np.zeros(3) + 1E6
    max_coords = np.zeros(3) - 1E6
    pdb_in = open(pdb_file, 'r')
    lines = pdb_in.readlines()
    pdb_in.close()
    for line in lines:
        data = line.split()
        if data[0] == "ATOM" or data[0] == "HETATM":
            name = data[2][0]
            if heavy and name == "H":
                continue  # Ignore hydrogens...
            atom_coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if type(angles) != type(None):
                atom_coords = rotate_coords(atom_coords, angles)
            # Pad coordinates
            for i in range(3):
                if atom_coords[i] > max_coords[i]:
                    max_coords[i] = atom_coords[i]
                if atom_coords[i] < min_coords[i]:
                    min_coords[i] = atom_coords[i]
    lengths = max_coords - min_coords
    volume = lengths[0] * lengths[1] * lengths[2]
    return volume


def write_pdb(input, output, angles, forward=True):
    """
    Write a PDB with rotated atoms into a new file
    forward variable gives the option of reversing a rotation...
    """
    pdb_in = open(input, 'r')
    in_lines = pdb_in.readlines()
    pdb_in.close()
    out_lines = []
    for line in in_lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            new_coords = rotate_coords(coords, angles, forward)
            new_line = line[:30] + '{:8.3f}{:8.3f}{:8.3f}'.format(new_coords[0], new_coords[1], new_coords[2])
            if len(line) > 54:
                new_line += line[54:]
            else:
                new_line += '\n'
        else:
            new_line = line
        out_lines.append(new_line)
    pdb_out = open(output, "w")
    for line in out_lines:
        pdb_out.write(line)
    pdb_out.close()
    degrees = angles * 180 / np.pi
    print("Structure from %s being written to %s after rotations of %.3f, %.3f and %.3f degrees about the x, y and z axes..."%(input, output, degrees[0], degrees[1], degrees[2]))


if __name__ == "__main__":
    # Haven't seen a need to write this as a stand-alone script yet - will likely just be used with the make_gcmcbox.py script
    print '\nDoes not currently operate as a stand-alone script\n'
    quit()

