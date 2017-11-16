"""
Script for network-based clustering of water arrangments
Intended as an alternative to the original clustering procedure employed to post-process GCMC simulations

This code is still highly experimental - the details have not been fully worked out...

Marley Samways
"""

import logging
import math

import numpy as np
from scipy.cluster import hierarchy

import simulationobjects

import copy
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from calc_clusters import _GetMolTemplate,_printpdb

def calc_distance(coords1, coords2):
    """
    Calculate distance between two sets of coordinates
    """
    vector = coords1 - coords2
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
    distance_matrix = np.zeros((len(coordlist1), len(coordlist2)))
    for i in range(len(coordlist1)):
        for j in range(len(coordlist2)):
            distance_matrix[i,j] = calc_distance(coordlist1[i], coordlist2[j])
    while len(pair_list) < len(coordlist1) and len(pair_list) < len(coordlist2):
        min_dist = 1E6
        for i in range(len(coordlist1)):
            already_counted = False
            for pair in pair_list:
                if pair[0] == i:
                    already_counted = True
            if already_counted: continue
            for j in range(len(coordlist2)):
                already_counted = False
                for pair in pair_list:
                    if pair[1] == j:
                        already_counted = True
                if already_counted: continue
                if distance_matrix[i,j] < min_dist:
                    min_dist = distance_matrix[i,j]
                    min_point = [i, j]
            pair_list.append(min_point)
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
    Find overlap of two lists..
    """
    overlap = []
    for i in list1:
        if i in list2:
            overlap.append(i)
    return overlap


def is_contained(list1, list2):
    """
    Check if list1 is contained within list2
    """
    if len(list1) > len(list2):
        return False
    elif len(list_overlap(list1, list2)) == len(list1):
        return True
    else:
        return False


def get_args():
    import argparse
    parser = argparse.ArgumentParser('Network-based clustering of hydration sites')
    parser.add_argument('-i', '--input', help='PDB file containing input frames. Currently accepts only one file', default='all.pdb')
    parser.add_argument('-m', '--molecule', help='Residue name of water molecules', default='WA1')
    parser.add_argument('-a', '--atom', help='Name of atom to take as molecule coordinates', default='O00')
    parser.add_argument('-s', '--skip', type=int, help='Number of frames to skip', default=0)
    parser.add_argument('-d', '--distance', help='Distance measure used: rmsd or max', default='rmsd')
    parser.add_argument('-l', '--linkage', help='Linkage method for hierarchical clustering', default='average')
    parser.add_argument('-c', '--cutoff', type=float, help='Distance cutoff to use during clustering', default=1.0)
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
    Nstar = np.mean(num_wat_list) 

    # STEP 2: Cluster frames within each bin, using RMSD as a measure of distance

    all_clusters = []
    
    for i in range(max_wats+1):
        if len(bins[i]) == 0: continue
        if len(bins[i]) == 1:
            all_clusters.append(bins[i])  # Cannot perform clustering if there is only one frame
	    continue
        dist_list = []
        for j in range(len(bins[i])):
            for k in range(j+1, len(bins[i])):
                if args.distance == 'rmsd':
                    dist_list.append(calc_rmsd(pdbfiles.pdbs[bins[i][j]].residues, pdbfiles.pdbs[bins[i][k]].residues))
                elif args.distance == 'max':
                    dist_list.append(calc_max_dist(pdbfiles.pdbs[bins[i][j]].residues, pdbfiles.pdbs[bins[i][k]].residues))
        # Perform hierarchical clustering
        tree = hierarchy.linkage(dist_list, method=args.linkage)
        # Use a distance cutoff to cut the tree
        clust_ids = hierarchy.fcluster(tree, t=args.cutoff, criterion='distance')
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
    principal_networks = []
    for i in range(num_clusters):
        cluster = all_clusters[i]
        min_dist = 1E6
        for j in range(len(cluster)):
            dist_list = []
            for k in range(len(cluster)):
                if args.distance == 'rmsd':
                    dist_list.append(calc_rmsd(pdbfiles.pdbs[cluster[j]].residues, pdbfiles.pdbs[cluster[k]].residues))
                elif args.distance == 'max':
                    dist_list.append(calc_max_dist(pdbfiles.pdbs[cluster[j]].residues, pdbfiles.pdbs[cluster[k]].residues))
            mean_dist = sum(dist_list) / len(dist_list)
            if mean_dist < min_dist:
                min_dist = mean_dist
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
        principal_networks.append(pdbfiles.pdbs[cluster[best_j]])
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
    # Not sure how best to do this yet...
    
    wat_coords = []
    wat_net_ids = []
    for i in range(len(principal_networks)):
        for j, wat in principal_networks[i].residues.iteritems():
            for atom in wat.atoms:
                if atom.name == "O00":
                    wat_coords.append(atom.coords)
                    wat_net_ids.append(i)
    wat_dist_list = []
    for i in range(len(wat_coords)):
        for j in range(i+1, len(wat_coords)):
            if wat_net_ids[i] == wat_net_ids[j]:
                wat_dist_list.append(1E6)
            else:
                wat_dist_list.append(calc_distance(wat_coords[i], wat_coords[j]))
    wat_tree = hierarchy.linkage(wat_dist_list, method='average')
    wat_clust_ids = hierarchy.fcluster(wat_tree, t=2.0, criterion='distance')
    wat_clusts = [[] for i in range(max(wat_clust_ids))]
    for i in range(len(wat_clust_ids)):
        clust_no = wat_clust_ids[i]-1
        wat_clusts[clust_no].append(wat_net_ids[i])
    pn_clusts = [[] for i in range(len(principal_networks))]
    for i in range(len(principal_networks)):
        for j in range(len(wat_clusts)):
            if i in wat_clusts[j]:
                pn_clusts[i].append(j)
    for i in range(len(principal_networks)):
        print("PN {} contains clusters\n\t{}".format(i+1, pn_clusts[i]))
    max_occupancy = max([sum([occupancies[j] for j in wat_clusts[i]]) for i in range(len(wat_clusts))])
    print("Maximum occupancy: {} %".format(max_occupancy))
    overlaps = [] # Lists of clusters which represent overlaps between PNs
    clust_sizes = []  # Numbers of clusters in each overlap
    sub_occupancies = [] # Occupancy of each subnetwork
    for i in range(len(pn_clusts)):
        for j in range(i+1, len(pn_clusts)):
            #print("Overlap between PNs {} and {}".format(i+1, j+1))
            overl = list_overlap(pn_clusts[i], pn_clusts[j])
            if not overl in overlaps:
                overlaps.append(overl)
                in_pns = []
                for k in range(len(pn_clusts)):
                    if is_contained(overl, pn_clusts[k]):
                        in_pns.append(k)
                clust_sizes.append(len(overl))
                sub_occupancies.append(sum([occupancies[k] for k in in_pns]))
    tested_until = 0  # Used to keep track of which sub-networks have been overlapped with PNs
    while max(sub_occupancies) < max_occupancy:
        # Keep finding overlaps until the highest occupancy unit has been found
        for i in range(tested_until, len(overlaps)):
            for j in range(len(pn_clusts)):
                if not is_contained(overlaps[i], pn_clusts[j]):
                    # Test the overlap against PNs which don't contain this entire unit
                    overl = list_overlap(overlaps[i], pn_clusts[j])
                    if not overl in overlaps:
                        overlaps.append(overl)
                        in_pns = [] 
                        for k in range(len(pn_clusts)):
                            if is_contained(overl, pn_clusts[k]):
                                in_pns.append(k)
                        clust_sizes.append(len(overl))
                        sub_occupancies.append(sum([occupancies[k] for k in in_pns]))
            tested_until += 1
    print("\n\nSub-networks from PN overlaps:")

    for i in range(len(overlaps)):
        print 'Subnetwork ',i,':'
        print("\t{}".format(overlaps[i]))
        in_pns = []
        for j in range(len(pn_clusts)):
            if is_contained(overlaps[i], pn_clusts[j]):
                in_pns.append(j)
        print("\t\tContains {} waters".format(len(overlaps[i])))
        clust_sizes.append(len(overlaps[i]))
        print("\t\tContained in PNs {}".format([j+1 for j in in_pns]))
        print("\t\tPresent in {:.2f} % of frames\n".format(sum([occupancies[j] for j in in_pns])))
        sub_occupancies.append(sum([occupancies[j] for j in in_pns]))
    # Add full PNs
    for i in range(len(pn_clusts)):
        if pn_clusts[i] in overlaps:
            continue
        overlaps.append(pn_clusts[i])
        clust_sizes.append(len(pn_clusts[i]))
        sub_occupancies.append(occupancies[i])

##### tree diagram
    print 'Plotting tree diagram from PNs and SubNs'
    zippeda = zip(range(1,len(overlaps)+1),sub_occupancies,clust_sizes,overlaps)
    zippeda.sort(key = lambda t: t[2])
    zipped = filter(lambda (a,b,c,d): b > 25.,zippeda)

    G=nx.Graph()
    sort_occ = []

   
    for row in zipped:
        sort_occ.append(row[1])

    additional_occ = [0.]*len(overlaps) 

    labels = {} 
    labels_back = {} 
    names = {}

    clust_sizes_filter = [c for a,b,c,d in zipped]
    occurence = {}
    clust_unique = list(set(clust_sizes_filter)) 
    for clust_size in clust_unique:
      occurence[clust_size]=clust_sizes_filter.count(clust_size)

    pos = {}
    previous = 0
    count = 0
    for i,a in enumerate(zipped):
      clust_size = zipped[i][2]
      if occurence[clust_size] == 1:
        pos[i+1]= (0,zipped[i][2])
      else:
        if clust_size != previous:  # if it is the first
          spacing = np.linspace(-7,7,occurence[clust_size]) 
          count = 0
        pos[i+1]= (spacing[count],zipped[i][2])
        previous = clust_size
        count+=1
      labels[i+1] = int(a[0])
      names[i+1] = str(' '.join( repr(e) for e in a[3]))
    G.add_nodes_from(pos)    
    labels_back = {int(v): k for k, v in labels.items()}

    edges = []
    for a,b in itertools.combinations(zipped,2):
      if all(x in b[3] for x in a[3]): 
        if a[2] == b[2]-1 : 
          edges.append((labels_back[a[0]],labels_back[b[0]]))
   
    nx.draw_networkx_nodes(G,pos,node_color=sort_occ,vmin = 0, vmax = 100, linewidths=0.0, cmap = plt.cm.Blues)
    nx.draw_networkx_edges(G,pos,edgelist=edges, edge_color='grey',width=3,alpha=0.5)
    nx.draw_networkx_labels(G,pos,labels,font_size=10)
   
    axes = plt.gca()
    axes.get_xaxis().set_visible(False)
    axes.set_yticks(np.arange(min(clust_sizes),max(clust_sizes)+1,1))
    
    spacing = 0.2
#    for each in names:
#      axes.text(0, min(clust_sizes)-0.2-spacing, str(each)+':  '+names[each],fontsize=10,horizontalalignment='center',color='green')
#      spacing+=0.2
    axes.set_ylim(min(clust_sizes)-spacing-0.2,max(clust_sizes)+1)
    axes.set_xlim(-8,8)
    #plt.plot(-8,Nstar,marker='*',color='orange',markersize=20)
    plt.plot([-8,8],[Nstar,Nstar],'--',alpha=0.5,color='orange')
    plt.savefig('tree.png')
    plt.gcf().clear()
    ### end of tree diagram --- take this directly from marleys orig....


# identifying core water molecules from largest occupancy sub_network
    num_of_core = 0.
    if max(sub_occupancies) >= 90.:
      print 'CORE WATER MOLECULES'
      print 'occupancy:', max(sub_occupancies)
      core = set(overlaps[sub_occupancies.index(max(sub_occupancies))]) 
      print 'water molecules:', ', '.join(str(s) for s in core) 
      num_of_core = len(core)
    else:
      print 'no persistent water molecules found. Site may be very disordered'
      core = set([])
   
# identifying non-core water molecules. This looks at the waters in the subnetworks with N* occupancies
    non_core_subnet=[] 
    for clust,waters in zip(clust_sizes,overlaps):
      if clust==math.floor(Nstar) or clust==math.ceil(Nstar):
        non_core_subnet.append(waters)
    non_core = [a for sublist in non_core_subnet for a in sublist]
    non_core = set(non_core)-set(core)
    print 'NON-CORE WATER MOLECULES'
    print 'water molecules:', ', '.join(str(s) for s in non_core) 

# this now finds the centroids of the clusters
    centroids = [[] for i in range(max(wat_clust_ids))] #centroids to be stored
    contained = [[] for i in range(max(wat_clust_ids))] #determine waters in each cluster
    for coords, ids in zip(wat_coords,wat_clust_ids):
      contained[ids-1].append(coords)
    
    for i,each in enumerate(contained): # determining centroid by average of all waters inside it
      centroids[i].append(np.mean(each,axis=0))


# rereading in data (re-write to avoid)
# below is taken from sub-network script to read in locations of ALL waters, not just PN ones
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    molpdb = _GetMolTemplate(args.input, args.molecule)
    num_frames = len(pdbfiles.pdbs)

    if len(pdbfiles.pdbs[0].residues) > 20.:
      quit()
    wat_list = []  # Store water molecules
    wat_whole = []  # Store water molecules
    frame_wat_ids = [[] for i in range(num_frames)]  # Store list ids of waters for each frame
    for i in range(num_frames):
        for j, wat in pdbfiles.pdbs[i].residues.iteritems():
            for atom in wat.atoms:
                if atom.name.upper() == 'O00':
                   wat_list.append(atom)
                   wat_whole.append(wat)
            frame_wat_ids[i].append(len(wat_list)-1)
 
    frame_clust_ids =copy.deepcopy(frame_wat_ids)  # converting water_ids to cluster ids in shape of frames ONLY NON CORE
    frame_clust_ids_all =copy.deepcopy(frame_wat_ids)  # converting water_ids to cluster ids in shape of frames CORE AND NONCORE
    frame_clust_dist =copy.deepcopy(frame_wat_ids)  # converting water_ids to cluster ids in shape of frames DISTANCE FROM EACH CENTROID
 
    for i,frame in enumerate(frame_wat_ids):
      for j,water in enumerate(frame):
        dist = [np.linalg.norm(x-wat_list[water].coords) for x in centroids]
        frame_clust_ids_all[i][j] = np.argmin(dist)
        frame_clust_dist[i][j] = np.min(dist)
        if np.argmin(dist) in non_core:
          frame_clust_ids[i][j] = np.argmin(dist)
        else:
          frame_clust_ids[i][j] = None


    flat_frame_clust_ids_all = [item for sublist in frame_clust_ids_all for item in sublist]
    all_occ = [0.]*len(centroids)
    for i,j in enumerate(range(len(centroids))):
      temp_occ = (flat_frame_clust_ids_all.count(j)/float(num_frames))*100
      if temp_occ > 100.:
        temp_occ = 100.
      all_occ[i] = temp_occ #calculating the overall occupancy of each cluster
      # need to write check that all_occ is never more than 100% occupied
      # a clustid should not be repeated in frame_clust_ids_all...



    for sub_id,clust in enumerate(non_core_subnet): # for each subnet of size N*
      best=[1E6,1E6]
      for i,(dist,frame) in enumerate(zip(frame_clust_dist,frame_clust_ids_all)): ### finding our best snapshot
        if all(x in frame for x in clust):
          distclust=0.
          for water in clust:
            distclust+= dist[frame.index(water)]
            if distclust <= best[1]:   #finding the frame with smallest overall distance (equivalent to RMSD)
              best[0] = i
              best[1] = distclust
      rep_struct_ids = []
      for water in clust:
        rep_struct_ids.append( frame_wat_ids[best[0]][frame_clust_ids_all[best[0]].index(water)])
      outfile = open('subnet'+str(sub_id)+'.pdb',"w") 
      for clustid,IDs in zip(clust,rep_struct_ids):
        coordinates = []
        for atom in wat_whole[IDs].atoms:
          coordinates.append( atom.coords)
        occ = float(all_occ[clustid])   # again, need to check that occ is not more than 100%
        _printpdb(np.asarray(coordinates),molpdb,clustid,occ,outfile)


    print 'Plotting correlation between non-core water molecules'

    G = nx.cycle_graph(len(non_core))
    pos = nx.circular_layout(G)   # arranges 'nodes' in a circle
    fig, ax = plt.subplots(1, 1, num=1)


    flat_frame_clust_ids = [item for sublist in frame_clust_ids for item in sublist]
    non_core_occ = [0.]*len(non_core)
    for i,j in enumerate(non_core):
      non_core_occ[i] = flat_frame_clust_ids.count(j)



    all_pairs = list(itertools.combinations(non_core, 2))
    observed_pairings=[]
    for frame in frame_clust_ids:
      for (i,j) in list(itertools.combinations(frame, 2)):
        observed_pairings.append((i,j))
    
    pair_occupancy = [0.]*len(all_pairs)
    weight = [0.]*len(all_pairs)
    for a,(i,j) in enumerate(all_pairs):
      pair_occupancy[a] = (observed_pairings.count((i,j)) + observed_pairings.count((j,i)))/float(num_frames)
      spontaneous_occupancy = (non_core_occ[list(non_core).index(i)]/float(num_frames))*(non_core_occ[list(non_core).index(j)]/float(num_frames))
      temp = 10*(pair_occupancy[a]-spontaneous_occupancy)
      if temp >=.5:
        weight[a] = temp 
      else:
        weight[a]= 0.
    for a,(i,j) in enumerate(all_pairs):
      all_pairs[a]= (list(non_core).index(i),list(non_core).index(j))

    ids = [x for x in list(non_core)]  # dictionary of included nodes
    labels = {}
    for i,j in enumerate(ids):
      labels[i]= str(j)

    nx.draw_networkx_nodes(G,pos,node_color=non_core_occ,cmap=plt.cm.Blues,vmin = 0., vmax = num_frames)
    nx.draw_networkx_edges(G,pos,edgelist=all_pairs,width=weight,edge_cmap=plt.cm.PRGn,alpha=0.5,edge_color=weight)
    nx.draw_networkx_labels(G,pos,labels,font_size=10)
    ax.set_axis_off()
    plt.savefig('correlations.png')
