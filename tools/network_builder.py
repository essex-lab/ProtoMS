"""
Test code for building a representative water network from GCMC simulation output

This code is very experimental!

Marley Samways
"""


import numpy as np
import argparse
from scipy.cluster import hierarchy
from scipy.stats import norm

import simulationobjects

def phi(a,b,frame_clust_ids):
    n11 = 0.0
    n01 = 0.0
    n10 = 0.0
    n00 = 0.0
    for frame in frame_clust_ids:
      if a in frame and b in frame:
        n11 += 1
      elif a in frame and b not in frame:
        n10 += 1
      elif a not in frame and b in frame:
        n01 += 1
      else:
        n00 +=1
    return (n11*n00 - n01*n10) / ((n11+n01)*(n11+n10)*(n00+n01)*(n00+n10))**0.5 

def get_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser('network_builder.py')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-d', '--distance', type=float, help='Distance cutoff for clustering. Default=3.0 Angs', default=3.0)
    parser.add_argument('-l', '--linkage', help='Linkage method for hierarchical clustering', default='average')
    parser.add_argument('-c', '--cutoff', type=float, help='Occupancy cutoffs for waters to be considered. Between 0 and 1.', default=0.0)
    parser.add_argument('-o', '--output', help='Stem for the PDB output', default='network')
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
    n_clusts = max(clust_ids)

    # count occupancy of each cluster
    clust_occ = [[i,list(clust_ids).count(i)] for i in range(1,n_clusts+1)] 
    clust_occ = sorted(clust_occ, key=lambda x:-x[1])
    # renumber the clusters based on this
    old_order = [x[0] for x in clust_occ]
    # this renumbers the cluster ids so that 1 is the most occupied...
    clust_ids_sorted = np.asarray([old_order.index(x)+1 for x in clust_ids])

    print("\n{} clusters identified:".format(n_clusts))
    clust_wat_ids = {}  # Store water IDs for each cluster
    for i in range(1, n_clusts+1):
        clust_wat_ids[i] = []
    for i in range(total_wats):
        clust = clust_ids_sorted[i]
        clust_wat_ids[clust].append(i)


# this prints out the average of each cluster 
    with open('clusts.pdb','w') as f: 
      for i in range(1, n_clusts+1):
          print("\tCluster {:2d}:\t{:3d}/{:3d} frames".format(i, len(clust_wat_ids[i]), n_frames))
      print("")
      for n in range(1,n_clusts+1):
          av_coords = [0.,0.,0.]
          count  = 0
          for i in clust_wat_ids[n]:
              for j in range(3):
                  av_coords[j] += coord_list[i][j]
              count += 1
          av_coords[0] /= count
          av_coords[1] /= count
          av_coords[2] /= count
          #print 'Cluster ', n, 'avg: ', av_coords
          f.write('HETATM{0:>5} {1:<4} {2:<4} {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.2f}         {9:>2}  \n'.format(n, 'O00', 'WA1',n, av_coords[0], av_coords[1], av_coords[2], clust_occ[n-1][1], clust_occ[n-1][1], 'O'))



    # Sort clusters
    clusts_sorted = []
    clust_occs = [len(clust_wat_ids[i]) for i in range(1, n_clusts+1)]
    while len(clusts_sorted) < n_clusts:
        max_occ = max([len(clust_wat_ids[i]) for i in range(1, n_clusts+1) if i not in clusts_sorted])
        for i in range(1, n_clusts+1):
            if i not in clusts_sorted:
                if np.isclose(clust_occs[i-1], max_occ):
                    clusts_sorted.append(i)

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
    # frame_clust_ids all of the simulation frames and the clusters observed in it


    # this finds the waters for which to do correlation analysis. Those that are on less than 100% of the time, and more than 25% of the time
    pos = []
    neg = []
    correlation_water_set = [(i,float(n)/float(n_frames)) for i,n in zip(clust_wat_ids,clust_occs) if n < n_frames and n >=0.25*n_frames]
    import itertools 
    # starting with pairwise, maybe possibility for n-body...
    for x,y in itertools.combinations(correlation_water_set,2):
      phi_coef = phi(x[0],y[0],frame_clust_ids)
      if phi_coef > 0.5:
        print
        print 'positive correlation'
        print 'between:',x[0],y[0]
        print 'phi:',np.round(phi_coef,2) 
        pos.append([x[0],y[0]])
      elif phi_coef < -0.5:
        print
        print 'negative correlation'
        print 'between:',x[0],y[0]
        print 'phi:',np.round(phi_coef,2) 
        neg.append([x[0],y[0]])
    with open('test.py','r') as i, open('good.py','w') as o:
      for line in i:
        if line[0:3] == 'pos':
          o.write('pos = '+str(pos)+'\n\n') 
        elif line[0:3] == 'neg':
          o.write('neg = '+str(neg)+'\n\n')
        else:
          o.write(line)
#    longlist = []
#    for pair in neg:
#      for y in pair:
#        longlist.append(y)
#    print set(longlist)
#    print all_in_neg
#    quit()

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
    
    # Need to identify which frames contain each network
    network_frame_ids = [[] for net in rep_networks]
    for i, network in enumerate(rep_networks):
        for j, frame in enumerate(frame_clust_ids):
            # Check that the network and frame contain the same clusters
            if len(network) == len(frame) and np.all([wat in frame for wat in network]):
                network_frame_ids[i].append(j)
    
    # Now need to compute cluster centres for each cluster
    clust_centres = {}
    for i in range(1, n_clusts+1):
        centre = np.zeros(3)
        for wat in clust_wat_ids[i]:
            centre += coord_list[wat]
        centre /= clust_occs[i-1]
        clust_centres[i] = centre
    
    # Now want to find a representative frame for each network
    network_rep_frames = []
    for i, network in enumerate(rep_networks):
        min_rmsd = 1E6  # Use smallest RMSD to identify the best frame
        for j in network_frame_ids[i]:
            frame = frame_wat_ids[j]
            sq_dists = []  # List of square differences in water positions
            for wat in frame:
                # Check which cluster each water in the frame corresponds to
                for k in range(1, n_clusts+1):
                    if wat in clust_wat_ids[k]:
                        sq_dists.append(np.square(np.linalg.norm(clust_centres[k]-coord_list[wat])))
                    else: continue
            rmsd = np.sqrt(sum(sq_dists))
            if rmsd < min_rmsd:
                min_rmsd = rmsd
                best_frame = j
        network_rep_frames.append(best_frame)
    
    # Print out networks
    print("{} representative network(s) built.".format(len(rep_networks)))
   
    # Write out representative frames to PDB files..
    print("Writing PDB output...")
    import time
    starttime = time.time() 
    for i, j in enumerate(network_rep_frames):
        # Write out water only for rep. frames
        pdbfiles.pdbs[j].write("{}-{}{}-wat.pdb".format(args.output, "0"*(2-len(str(i+1))), i+1))
        # Read in whole system - only want to read a single frame to save time/memory
        pdbfiles2 = simulationobjects.PDBSet()
        pdbfiles2.read(args.input, skip=args.skip+j, readmax=1)
        pdbfiles2.pdbs[0].write("{}-{}{}-all.pdb".format(args.output, "0"*(2-len(str(i+1))), i+1))
    print("{} seconds to read in all PDB data and then write out frames".format(time.time()-starttime))
        
    print("Done!\n")

