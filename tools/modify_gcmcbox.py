"""
Script to modify GCMC box during system setup
Requires that a box has already been created
"""

import numpy as np

import simulationobjects
from make_gcmcbox import print_bequil


def translate_box(origin, translation):
    """
    Function to shift the box through space by some vector
    """
    new_origin = np.zeros(3)
    for i in range(3):
        new_origin[i] = origin[i] + translation[i]
    return new_origin


def pad_box(lengths, paddings):
    """
    Function to pad the x, y and z-axes by some values input by the user
    """
    new_lengths = np.zeros(3)
    for i in range(3):
        new_lengths[i] = lengths[i] + 2*paddings[i]
    return new_lengths


def read_box(boxpdb):
    """
    Read origin and dimensions of a GCMC box...
    """
    input_file = open(args.input, 'r')
    input_lines = input_file.readlines()
    input_file.close()

    for line in input_lines[:3]:
        data = line.split()
        if data[1] == "CENTER":
            origin = np.array([float(data[5]), float(data[6]), float(data[7])])
        elif data[1] == "DIMENSIONS":
            lengths = np.array([float(data[5]), float(data[6]), float(data[7])])
    return origin, lengths



def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="Script to translate and pad along axes of a GCMC box")
    parser.add_argument('-i', '--input', help="the name of the PDB-file containing the original GCMC box", default='gcmc_box.pdb')
    parser.add_argument('-o', '--output', help="the name of the PDB-file containing the modified GCMC box", default='gcmc_box2.pdb')
    parser.add_argument('-tx', '--xtrans', type=float, help="Translation along x-axis", default=0.0)
    parser.add_argument('-ty', '--ytrans', type=float, help="Translation along y-axis", default=0.0)
    parser.add_argument('-tz', '--ztrans', type=float, help="Translation along z-axis", default=0.0)
    parser.add_argument('-px', '--xpad', type=float, help='Extra padding along the x-axis', default=0.0)
    parser.add_argument('-py', '--ypad', type=float, help='Extra padding along the y-axis', default=0.0)
    parser.add_argument('-pz', '--zpad', type=float, help='Extra padding along the z-axis', default=0.0)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    
    # Read in box parameters
    origin_in, lengths_in = read_box(args.input)
    
    # Translate and pad the box
    origin_out = translate_box(origin_in, [args.xtrans, args.ytrans, args.ztrans])
    lengths_out = pad_box(lengths_in, [args.xpad, args.ypad, args.zpad])
    
    # Calculate volume and print Bequil value
    volume = lengths_out[0] * lengths_out[1] * lengths_out[2]
    print_bequil(volume)

    # Write new box to a file
    box_out = {"center":origin_out, "len":lengths_out}
    simulationobjects.write_box(args.output, box_out)






