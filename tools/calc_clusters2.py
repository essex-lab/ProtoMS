"""
calc_clusters2.py
Script to cluster water observations from a GCMC simulation.
Intended to replace calc_clusters.py

Marley Samways
Hannah Bruce Macdonald
Chris Cave-Ayland
Gregory Ross
"""

import argparse
import numpy as np
import six
from scipy.cluster import hierarchy

from protomslib import simulationobjects


def get_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser('calc_clusters2.py')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-c', '--cutoff', type=float, help='Distance cutoff for clustering. Default=3.0 Angs', default=3.0)
    parser.add_argument('-l', '--linkage', help='Linkage method for hierarchical clustering', default='average')
    parser.add_argument('-o', '--output', help='Filename for the PDB output', default='clusts.pdb')
    args = parser.parse_args()
    return args


def get_distances(frame_wat_ids, coord_list):
    """
    Calculate the distances between all water observations as a flatted 1D distance matrix

    Parameters
    ----------
    frame_wat_ids : list
        List of lists containing the water IDs present in each frame
    coord_list : list
        List of oxygen coordinates for each water molecule observed

    Returns
    -------
    dist_list : list
        List of distances between water molecules as a 1D matrix (flattened from 2D)
    """
    print("Calculating water-water distances...")
    n_wats = len(coord_list)
    dist_list = []
    for i in range(n_wats):
        frame_last = True  # Check if last water was in the same frame
        for j in range(i+1, n_wats):
            same_frame = False
            if frame_last:
                for frame in frame_wat_ids:
                    if i in frame and j in frame:
                        same_frame = True
                        break
            if same_frame:
                dist_list.append(1E8)
            else:
                dist_list.append(np.linalg.norm(coord_list[i]-coord_list[j]))
                frame_last = False  # Stops checking if subsequent j values are in the same frame
    return dist_list


def sort_clusters(clust_ids):
    """
    Calculate cluster occupancies and order them based on these values

    Parameters
    ----------
    clust_ids : list
        List indicating which cluster each observation belongs to

    Returns
    -------
    clust_ids_sorted : list
        Cluster IDs sorted by occupancy
    clust_occ : list
        List of occupancies for each cluster
    """
    n_clusts = max(clust_ids)
    # Count occupancy of each cluster
    clust_occ = [[i, list(clust_ids).count(i)] for i in range(1, n_clusts+1)]
    clust_occ = sorted(clust_occ, key=lambda x:-x[1])
    # Renumber the clusters based on this
    old_order = [x[0] for x in clust_occ]
    clust_occs = [x[1] for x in clust_occ]
    # This renumbers the cluster ids for each water so that 1 is the most occupied...
    clust_ids_sorted = np.asarray([old_order.index(x)+1 for x in clust_ids])
    return clust_ids_sorted, clust_occs


def write_clusts(filename, clust_wat_ids, n_clusts, n_frames, clust_occs, all_wats):
    """
    Write the average oxygen positions of each cluster to a PDB file

    Parameters
    ----------
    filename : str
        Name of PDB output file
    clust_wat_ids : dict
        Dictionary indicating which water observations are present in each cluster
    n_clusts : int
        Total number of clusters
    n_frames : int
        Total number of frames
    clust_occs : list
        List containing occupancies of each cluster
    all_wats : list
        List of all water observations, each an instance of 

    Returns
    -------
    centres : dict
        Dictionary containing the average coordinates of each cluster
    """
    # Print cluster occupancy data to screen
    for i in range(1, n_clusts+1):
        print("\tCluster {:2d}:\t{:3d}/{:3d} frames".format(i, len(clust_wat_ids[i]), n_frames))
    print("")
    # Calculate representative coordinates and write to PDB
    with open(filename, 'w') as f:
        f.write("REMARK Average cluster locations written by calc_clusters2.py\n")
        for n in range(1, n_clusts+1):
            # Calculate average coordinates
            av_coords = np.zeros(3)
            for i in clust_wat_ids[n]:
                av_coords += coord_list[i]
            av_coords /= len(clust_wat_ids[n])
            # Find the closest observation to the centre and use these coordinates as a representative
            min_dist = 1E6
            for i in clust_wat_ids[n]:
                if np.linalg.norm(coord_list[i]-av_coords) < min_dist:
                    rep_wat = all_wats[i]
            # Write to file
            for i, atom in enumerate(rep_wat.atoms):
                atom_id = (n - 1) * len(rep_wat.atoms) + i + 1
                f.write('ATOM  {0:>5} {1:<4} {2:<4} {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.2f}\n'.format(atom_id, atom.name, 'WA1', n,
                                                                                                                      atom.coords[0], atom.coords[1], atom.coords[2],
                                                                                                                      clust_occs[n-1]/float(n_frames), clust_occs[n-1]))
            f.write('TER\n')
        f.write('END')
    return None


if __name__ == "__main__":
    # Read command line arguments
    args = get_args()
    
    # Read PDB input data - only water molecules
    print('\nReading PDB data...')
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    n_frames = len(pdbfiles.pdbs)
    print('{} frames of PDB data read in.\n'.format(n_frames))
    
    # Store all water molecules and oxygen coordinates
    wat_list = []
    coord_list = []
    frame_wat_ids = [[] for i in range(n_frames)]
    for i in range(n_frames):
        for j, wat in six.iteritems(pdbfiles.pdbs[i].residues):
            wat_list.append(wat)
            for atom in wat.atoms:
                if atom.name == args.atom:
                    coord_list.append(atom.coords)
            frame_wat_ids[i].append(len(wat_list)-1)
    total_wats = len(wat_list)

    # Calculate a 1D distance matrix between all water observations
    dist_list = get_distances(frame_wat_ids, coord_list)
    
    # Cluster waters into discrete hydration sites
    print("Clustering water sites...")
    tree = hierarchy.linkage(dist_list, method=args.linkage)
    clust_ids = hierarchy.fcluster(tree, t=args.cutoff, criterion='distance')
    n_clusts = max(clust_ids)
    print("\n{} clusters identified:".format(n_clusts))

    # Renumber clusters based on occupancy
    clust_ids_sorted, clust_occs = sort_clusters(clust_ids)

    # Identify which water observations are present in each cluster
    clust_wat_ids = {}  # Store water IDs for each cluster
    for i in range(1, n_clusts+1):
        clust_wat_ids[i] = []
    for i in range(total_wats):
        clust = clust_ids_sorted[i]
        clust_wat_ids[clust].append(i)

    # Prints out the representative of each cluster to PDB
    write_clusts("{}".format(args.output), clust_wat_ids,
                 n_clusts, n_frames, clust_occs, wat_list)


