"""
Script to find long live water units in the hydration patterns from GCMC runs

Marley Samways
"""

import logging

import numpy as np
from scipy.cluster import hierarchy

import simulationobjects
from calc_clusters import _GetMolTemplate,_printpdb


def calc_distance(wat1, wat2):
    """
    Calculate distance between two sets of coordinates
    """
    vector = wat1.coords - wat2.coords
    dist_sq = sum(np.square(vector))
    return np.sqrt(dist_sq)


def find_pairs(watlist1, watlist2):
    """
    Match each water in a set to a water in another set
    """
    # Store water coordinates via representative atoms
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
    return pair_list, coordlist1, coordlist2


def calc_rmsd(watlist1, watlist2):
    """
    Calculate the RMSD of two sets of water molecules
    """
    # Match each the waters in each set to form pairs
    pair_list, coordlist1, coordlist2 = find_pairs(watlist1, watlist2)
    # Calculate RMSD
    sq_dist_list = []
    for pair in pair_list:
        i = pair[0]
        j = pair[1]
        sq_dist_list.append(np.square(calc_distance(coordlist1[i], coordlist2[j])))
    mean_sq_dist = sum(sq_dist_list) / len(sq_dist_list)
    rmsd = np.sqrt(mean_sq_dist)
    return rmsd


def calc_max_dist(watlist1, watlist2):
    """
    Find maximum distance between two waters for each of the pairs
    """
    # Match the waters in each set
    pair_list, coordlist1, coordlist2 = find_pairs(watlist1, watlist2)
    # Find maximum distance
    dist_list = []
    for pair in pair_list:
        i = pair[0]
        j = pair[1]
        dist_list.append(np.square(calc_distance(coordlist1[i], coordlist2[j])))
    max_dist = max(dist_list)
    return max_dist


def list_overlap(list1, list2):
    """
    Find overlap between two lists
    """
    if len(list1) > len(list2):
        long_list = list1
        short_list = list2
    else:
        long_list = list2
        short_list = list1
    overlap_list = []
    for i in long_list:
        if i in short_list:
            overlap_list.append(i)
    return overlap_list


def is_contained(list1, list2):
    """
    Checks if list1 is fully contained within list2
    """
    if len(list1) > len(list2):
        return False
    contained = True
    for x in list1:
        if not x in list2:
            contained = False
            break
    return contained

def get_representative_struct(net_id,frames,frame_wat_ids,watid_to_clustid):
  coordinates = [[] for i in range (len(net_id))]
  for frame in frames:
    for water in frame_wat_ids[frame]:
      if watid_to_clustid[water] in net_id:
	coordinates[net_id.index(watid_to_clustid[water])].append(wat_list[water].coords)
  coordinates = np.asarray(coordinates)
  average_network = np.mean(coordinates,axis=1)
  diff = 1E6
  for i,frame in enumerate(frames):
    snapshot = coordinates[:,i,:]
    if np.mean(np.abs(snapshot-average_network)) < diff:
      diff = np.mean(np.abs(snapshot-average_network))
      best_frame = frame 
  ids_in_best_frame = []
  for water in frame_wat_ids[best_frame]:
    if watid_to_clustid[water] in net_id:
      ids_in_best_frame.append(water)
  return ids_in_best_frame

def get_args():
    import argparse
    parser = argparse.ArgumentParser('Network-based clustering of hydration sites')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-l', '--linkage', help='Linkage method for hierarchical clustering', default='average')
    parser.add_argument('-c', '--cutoff', type=float, help='Distance cutoff to use during clustering', default=2.0)
    parser.add_argument('-occ', '--occupancy', type=float, help='Percentage of frames desired for network occupancies', default=40.0)
    parser.add_argument('-out', '--output', help='Stem for the output networks', default='subNet_')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    
    # Read PDB input data
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    molpdb = _GetMolTemplate(args.input, args.molecule)
    num_frames = len(pdbfiles.pdbs)
    print '\n{} frames of PDB data read in.\n'.format(num_frames)
    
    # STEP 1: Load all waters from all frames into a list, noting which frame each water came from

    wat_list = []  # Store water molecules
    wat_whole = []  # Store water molecules
    frame_wat_ids = [[] for i in range(num_frames)]  # Store list ids of waters for each frame
    for i in range(num_frames):
        for j, wat in pdbfiles.pdbs[i].residues.iteritems():
	    wat_whole.append(wat)
            for atom in wat.atoms:
                if atom.name.upper() == 'O00':
                    wat_list.append(atom)
            frame_wat_ids[i].append(len(wat_list)-1)

    # STEP 2: construct a 1D distance matrix of all water molecules.
    # If two waters are within the same frame, the distance is set to 1E6
    print 'Constructing distance matrix and checking which waters are in the same frame...'
    dist_list = []
    for i in range(len(wat_list)):
        frame_last = True  # Checks if the last water was in the same frame
        for j in range(i+1, len(wat_list)):
            same_frame = False
            if frame_last:
                for frame in frame_wat_ids:
                    if i in frame and j in frame:
                        same_frame = True
                        break
            if same_frame:
                dist = 1E6
            else:
                dist = calc_distance(wat_list[i], wat_list[j])
                frame_last = False  # Stops checking if subsequent waters are in the same frame for this i...
            dist_list.append(dist)

    # STEP 3: Perform hierarchical clustering to bin observations into hydration sites

    # Perform hierarchical clustering to get a linkage tree
    print 'Creating hierarchy...'
    tree = hierarchy.linkage(dist_list, method=args.linkage)
    # Use a distance cutoff to cut the tree
    print 'Cutting tree...\n'
    clust_ids = hierarchy.fcluster(tree, t=args.cutoff, criterion='distance')
#    k = 0 
#    for i, j in zip(clust_ids,wat_list):
#      print k
#      k += 1
#      print i
#      print j.coords
   
    cluster_wat_ids = [[] for i in range(max(clust_ids))]  # List to store ids of waters contained in each cluster
    print '{} hydration sites identified'.format(max(clust_ids))
    clust_nos = clust_ids - 1 
    for i in range(len(clust_ids)):
        clust_no = clust_ids[i] - 1
        cluster_wat_ids[clust_no].append(i)
    cluster_frame_ids = [[] for i in range(max(clust_ids))]
    for i in range(max(clust_ids)):
        for wat_id in cluster_wat_ids[i]:
            for j in range(len(frame_wat_ids)):
                if wat_id in frame_wat_ids[j]:
                    cluster_frame_ids[i].append(j)
                    break

    # STEP 4: Build sub-networks from hydration sites, provided that the sites are observed together a sufficient proporttion of the time

    network_ids = []  # Stores cluster ID list for each network formed
    network_frames = []
    for i in range(max(clust_ids)):
        if len(cluster_frame_ids[i])/float(num_frames) >= args.occupancy/100:
            network_ids.append([i])
            network_frames.append(cluster_frame_ids[i])
    merging = True  # Used to indicate that networks are being merged...
    while merging:
        merging = False
        for i in range(len(network_ids)):
            for j in range(i+1, len(network_ids)):
                if len(list_overlap(network_ids[i], network_ids[j])) != 0:
                    continue  # Makes sure that no waters are added to a network twice
                new_ids = []
                comb_ids = []
                for k in network_ids[i]:
                    comb_ids.append(k)
                for k in network_ids[j]:
                    comb_ids.append(k)
                for x in comb_ids:
                    if not x in new_ids:
                        new_ids.append(x)
                # Make sure the merged network is not contained elsewhere
                already_done = False
                for k in range(len(network_ids)):
                    if len(list_overlap(new_ids, network_ids[k])) == len(new_ids):
                        already_done = True
                        break
                if already_done:
                    continue
                new_frames = list_overlap(network_frames[i], network_frames[j])
                if len(new_frames)/float(num_frames) >= args.occupancy/100:
                    new_ids = []
                    for k in network_ids[i]:
                        new_ids.append(k)
                    for k in network_ids[j]:
                        new_ids.append(k)
                    network_ids.append(new_ids)
                    network_frames.append(new_frames)
            if merging:
                continue

    # STEP 5: Remove redundancy by filtering out sub-netowrks which are fully contained in larger sub-networks

    repeated_units = []  # Stores subnetwork IDs which are contained in other networks
    for i in range(len(network_ids)):
        for j in range(len(network_ids)):
            if i == j:
                continue
            if is_contained(network_ids[i], network_ids[j]):
                repeated_units.append(network_ids[i])
                break
    print '\n'
   # Print details of each network

    occupancy = []
    for frames in network_frames:
      occupancy.append(len(frames))
    collated = zip(network_ids,network_frames,occupancy)
    collated.sort(key = lambda t: t[2],reverse=True) # sorting networks from most to least populated

    watid_to_clustid = dict(zip( range(len(clust_nos)), clust_nos )) #-1 as clust_ids starts at 1 and network_ids start at 0

    j = 0
    for net_id, frames, occ in collated:
	if net_id in repeated_units:
	  continue
        j += 1
        print 'Sub-network {}:'.format(j) 
        print '\tContains clusters: {}'.format(net_id)
        print '\tPresent in {:.1f} % of frames'.format(100*occ/float(num_frames))
        print '\n'
	rep_struct_ids = get_representative_struct(net_id,frames,frame_wat_ids,watid_to_clustid)
  	outfile = open(args.output+str(j)+'.pdb',"w")
        for IDs in rep_struct_ids:
	  coordinates = []
	  for atom in wat_whole[IDs].atoms:
	    coordinates.append( atom.coords)
          _printpdb(np.asarray(coordinates),molpdb,watid_to_clustid[IDs],100.,outfile)

#    watid_to_clustid = dict(zip( range(len(clust_ids)), clust_ids ))
#    print watid_to_clustid
#    print conserved_frame_ids
#    for cluster in conserved_frame_ids: #looping over cluster
#      print cluster
#      for frame in cluster: # for each 
#	print 'frame ID', frame
#	print 'waters in frame',frame_wat_ids[frame]	
#	for water in frame_wat_ids[frame]:
#	  print water
#	  print watid_to_clustid[water]
