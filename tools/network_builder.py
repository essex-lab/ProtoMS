"""
Test code for building a representative water network from GCMC simulation output

This code is very experimental!

Marley Samways
"""


import numpy as np
import argparse
from scipy.cluster import hierarchy

import simulationobjects


def get_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser('network_builder.py')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-d', '--distance', type=float, help='Distance cutoff for clustering. Default=2.0 Angs', default=2.0)
    parser.add_argument('-l', '--linkage', help='Linkage method for hierarchical clustering', default='average')
    parser.add_argument('-c', '--cutoff', type=float, help='Occupancy cutoffs for waters to be considered. Between 0 and 1.', default=0.0)
    parser.add_argument('-o', '--output', help='Stem for the PDB output', default='output')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Read command line arguments
    args = get_args()
    
    # Read PDB input data
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    n_frames = len(pdbfiles.pdbs)
    print '\n{} frames of PDB data read in.\n'.format(n_frames)
    
    # Store all water molecules and oxygen coordinates
    wat_list = []
    coord_list = []
    frame_wat_ids = [[] for i in range(n_frames)]
    for i in range(n_frames):
        for j, wat in pdbfiles.pdbs[i].residues.iteritems():
            wat_list.append(wat)
            for atom in wat.atoms:
                if atom.name == args.atom:
                    coord_list.append(atom.coords)
            frame_wat_ids[i].append(len(wat_list)-1)
    total_wats = len(wat_list)
    
    # Calculate a one-dimensionalised distance matrix of all water observations
    print("Calculating distance matrix...")
    dist_list = []
    for i in range(total_wats):
        frame_last = True  # Check if last water was in the same frame
        for j in range(i+1, total_wats):
            same_frame = False
            if frame_last:
                for frame in frame_wat_ids:
                    if i in frame and j in frame:
                        same_frame = True
                        break
            if same_frame:
                dist_list.append(1E6)
            else:
                dist_list.append(np.linalg.norm(coord_list[i]-coord_list[j]))
                frame_last = False # Stops checking if subsequent j values are in the same frame
    
    # Cluster waters into discrete hydration sites
    print("Clustering...")
    tree = hierarchy.linkage(dist_list, method=args.linkage)
    clust_ids = hierarchy.fcluster(tree, t=args.distance, criterion='distance')
    
    clust_wat_ids = {}  # Store water IDs for each cluster
    n_clusts = max(clust_ids)
    for i in range(1, n_clusts+1):
        clust_wat_ids[i] = []
    for i in range(total_wats):
        clust = clust_ids[i]
        clust_wat_ids[clust].append(i)
    
    print("\nClusters:")
    for i in range(1, n_clusts+1):
        print("\tCluster {:2d}:\t{:3d}/{:3d} frames".format(i, len(clust_wat_ids[i]), n_frames))
    print("")
    
    # Sort clusters
    clusts_sorted = []
    clust_occs = [len(clust_wat_ids[i]) for i in range(1, n_clusts+1)]
    while len(clusts_sorted) < n_clusts:
        max_occ = max([len(clust_wat_ids[i]) for i in range(1, n_clusts+1) if i not in clusts_sorted])
        for i in range(1, n_clusts+1):
            if i not in clusts_sorted:
                if np.isclose(clust_occs[i-1], max_occ):
                    clusts_sorted.append(i)
    #print("Clusters in order of occupancy:")
    #print("\t{}\n".format(clusts_sorted))
    
    # Check which clusters are ever observed together
    clust_frame_ids = {}  # Store frame IDs for each cluster
    frame_clust_ids = [[] for i in range(n_frames)]
    for i in range(1, n_clusts+1):
        clust_frame_ids[i] = []
        for wat_id in clust_wat_ids[i]:
            for j, frame in enumerate(frame_wat_ids):
                if wat_id in frame:
                    clust_frame_ids[i].append(j)
                    frame_clust_ids[j].append(i)

    # Build a set of networks starting with the most occupied sites
    rep_networks = []
    left_out = [i for i in clusts_sorted if clust_occs[i-1] > args.cutoff * n_frames]
    while len(left_out) > 0:
        network = [] 
        # Start the network by including those which have been left out..
        for i in left_out + clusts_sorted:
            if i in network: continue  # Make sure not to include the same water twice...
            suitable = False  # Checks that all waters in the network have been observed together at least once
            # Check that all waters are observed in the same frame at least once
            for frame in frame_clust_ids:
                if i in frame and np.all([j in frame for j in network]):
                    suitable = True
            if suitable:
                network.append(i)
        # Add this network to the list
        rep_networks.append([i for i in clusts_sorted if i in network])  # Sort the cluster IDs in the network
        left_out = []
        # Check which occupied sites have not been yet included in a network.
        for i in clusts_sorted:
            if clust_occs[i-1] > args.cutoff * n_frames and np.all([i not in net for net in rep_networks]):
                left_out.append(i)
    
    # Need to identify which frames represent each network
    
    # Print out networks
    print("Representative network(s) built:")
    for network in rep_networks:
        print("\t{}\n".format(network))
    

