"""
Script for network-based clustering of water arrangments
Intended as an alternative to the original clustering procedure employed to post-process GCMC simulations

This code is still highly experimental - the details have not been fully worked out...

Marley Samways
"""

import logging

import numpy as np
from scipy.cluster import hierarchy

import simulationobjects


def calc_distance(coords1, coords2):
    """
    Calculate distance between two sets of coordinates
    """
    vector = coords1 - coords2
    dist_sq = sum(np.square(vector))
    return np.sqrt(dist_sq)


def calc_rmsd(watlist1, watlist2, atom):
    """
    Calculate the RMSD of two sets of water molecules
    """
    # Store coordinates of waters via representative atoms
    coordlist1 = []
    for i, wat in watlist1.iteritems():
        for atom in wat.atoms:
            if atom.name == "O00":
                coordlist1.append(atom.coords)
    coordlist2 = []
    for i, wat in watlist2.iteritems():
        for atom in wat.atoms:
            if atom.name == "O00":
                coordlist2.append(atom.coords)
    # Store a list of matched waters
    pair_list = []
    for i in range(len(coordlist1)):
        closest_dist = 1E6
        for j in range(len(coordlist2)):
            already_paired = False
            for pair in pair_list:
                if pair[1] == j:
                    already_paired = True
                    break
            if already_paired:
                continue
            dist = calc_distance(coordlist1[i], coordlist2[j])
            if dist < closest_dist:
                closest_dist = dist
                closest_j = j
        if closest_dist < 100:
            pair_list.append([i, closest_j])
    sq_dist_list = []
    for pair in pair_list:
        i = pair[0]
        j = pair[1]
        sq_dist_list.append(np.square(calc_distance(coordlist1[i], coordlist2[j])))
    mean_sq_dist = sum(sq_dist_list) / len(sq_dist_list)
    rmsd = np.sqrt(mean_sq_dist)
    return rmsd


def get_args():
    import argparse
    parser = argparse.ArgumentParser('Network-based clustering of hydration sites')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-rmsd', '--rmsdcut', type=float, help='RMSD cutoff to use during clustering', default=1.0)
    parser.add_argument('-o', '--output', help='Stem for the output networks', default='principalNetwork_')
    parser.add_argument('--plot', help='Plot the probability distributions of each of the principal networks.', action='store_true')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    
    # Read PDB input data
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    num_frames = len(pdbfiles.pdbs)
    print '\n{} frames of PDB data read in.\n'.format(num_frames)
    
    # STEP 1: Bin frames by number of waters, so that all frames containing the same number are lumped together

    # Find max number of waters in order to determine how many bins to create
    num_wat_list = []
    for i in range(num_frames):
        num_wat_list.append(len(pdbfiles.pdbs[i].residues))
    max_wats = max(num_wat_list)

    # For each bin, load the frame ids corresponding to the appropriate number of waters    
    bins = [[] for i in range(max_wats+1)]
    for i in range(num_frames):
        num_wats = len(pdbfiles.pdbs[i].residues)
        bins[num_wats].append(i)

    for i in range(max_wats+1):
        if len(bins[i]) == 0: continue
        print '{:>3} frames contain {} waters'.format(len(bins[i]), i)

    # STEP 2: Cluster frames within each bin, using RMSD as a measure of distance

    all_clusters = []
    
    for i in range(max_wats+1):
        if len(bins[i]) == 0: continue
        rmsd_list = []
        for j in range(len(bins[i])):
            for k in range(j+1, len(bins[i])):
                rmsd_list.append(calc_rmsd(pdbfiles.pdbs[bins[i][j]].residues, pdbfiles.pdbs[bins[i][k]].residues, args.atom))
        # Perform average-linkage hierarchical clustering
        tree = hierarchy.linkage(rmsd_list, method='average')
        # Use a distance cutoff to cut the tree
        clust_ids = hierarchy.fcluster(tree, t=args.rmsdcut, criterion='distance')
        num_new_clusts = max(clust_ids)
        new_clusts = [[] for j in range(num_new_clusts)]
        for j in range(len(bins[i])):
            clust_no = clust_ids[j] - 1
            new_clusts[clust_no].append(bins[i][j])
        for cluster in new_clusts:
            all_clusters.append(cluster)

    num_clusters = len(all_clusters)
    print '\n{} principal networks identified\n'.format(num_clusters)

    ordered_clusters = []  # Order clusters by occupancy

    # Puts most occupied networks in order
    while len(all_clusters) > 0:
        max_occupancy = 0
        for i in range(len(all_clusters)):
            occupancy = len(all_clusters[i])
            if occupancy > max_occupancy:
                max_occupancy = occupancy
                best_i = i
        ordered_clusters.append(all_clusters[best_i])
        all_clusters = all_clusters[:best_i] + all_clusters[best_i+1:]

    all_clusters = ordered_clusters

    # STEP 3: Use clusters to generate 'principal networks'

    occupancies = []

    for i in range(num_clusters):
        print 'Principal Network {}:'.format(i+1)
        print '    Number of waters:\t{}'.format(len(pdbfiles.pdbs[all_clusters[i][0]].residues))
        cluster_occupancy = len(all_clusters[i]) * 100.0 / num_frames
        occupancies.append(cluster_occupancy)
        print '    Frame occupancy: \t{:.2f} %\n'.format(cluster_occupancy)

    print 'Writing networks to PDB files...\n'
    # Write networks to PDB files - choose a representative frame as that which has the smallest RMSD with all others
    for i in range(num_clusters):
        cluster = all_clusters[i]
        min_rmsd = 1E6
        for j in range(len(cluster)):
            rmsd_list = []
            for k in range(len(cluster)):
                rmsd_list.append(calc_rmsd(pdbfiles.pdbs[cluster[j]].residues, pdbfiles.pdbs[cluster[k]].residues, args.atom))
            mean_rmsd = sum(rmsd_list) / len(rmsd_list)
            if mean_rmsd < min_rmsd:
                min_rmsd = mean_rmsd
                best_j = j
        num_wats = len(pdbfiles.pdbs[cluster[best_j]].residues)
        filename = args.output
        if num_wats < 10:
            filename += '0'
        filename += str(num_wats) + '_'
        if i + 1 < 100:
            filename += '0'
        if i + 1 < 10:
            filename += '0'
        filename += str(i+1) + '.pdb'
        pdbfiles.pdbs[cluster[best_j]].write(filename)

    if args.plot:
        # Plot the occupancies of the principal networks, in order
        import matplotlib.pyplot as plt
        ids = [i for i in range(num_clusters)]
        fig, ax = plt.subplots()
        bar_width = 1.0
        opacity = 1.0
        bar_pos = [ids[i]+0.5*bar_width for i in range(len(ids))] 
        bars = plt.bar(bar_pos, occupancies, bar_width, alpha=opacity, color='red')
        plt.xlabel('Principal network ID')
        plt.ylabel('Occupancy / %')
        plt.xlim(0.5, num_clusters+0.5)
        plt.ylim(0, 100)
        plt.show()
        

    # STEP 4: Condense the principal network data into a smaller number of more meaningful networks
    # Not sure how to do this yet...


