"""
Analysis code for building representative water network(s) from GCMC simulation output

This code is still under development!

Marley Samways
Hannah Bruce Macdonald
"""

import argparse
import itertools
import time
import numpy as np
from scipy.cluster import hierarchy

import simulationobjects


def separate_networks(check,phis,phi_cutoff):
    combs = []
    for i in range(1, len(check)):
        els = [list(x) for x in itertools.combinations(check,i)]
        for each in els:
            combination = [each, [x for x in check if x not in each]]
            if combination not in combs and [[x for x in check if x not in each],each] not in combs:
                combs.append(combination)
    
    bestscore = -100.
    bestcombs = []
    for x in combs:
        score = 0.
        for a in itertools.combinations(x[0],2):
            if phis[a] < -phi_cutoff:
                score -= 100 
            else:
                score+= phis[a]
        for a in itertools.combinations(x[1],2):
            if phis[a] < -phi_cutoff:
                score -= 100
            else:
                score+= phis[a]
        if score > bestscore:
            bestscore = score
            bestcombs = x
            
    if bestscore == -100.:
        print 'not satisfied with 2 networks'
    print 'this is the best network setup:',bestcombs
    print 'with this score',bestscore
    return bestcombs


def network_alt(check,phis):
    print check
    all_networks = []
    for a in check:
        this_network = [a]
	others = [x for x in check if x != a]
        for b in others:
            if a<b:
         	if phis[(a,b)] > -0.27:
    	            this_network.append(b)
	    else:
         	if phis[(b,a)] >= -0.27:
    	            this_network.append(b)
	all_networks.append(sorted(this_network))
    print all_networks
    print len(all_networks)
    print len(set(tuple(net) for net in all_networks))

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
    parser.add_argument('-p', '--phi', type=float, help='Magnitude of phi to indicate correlation. Between 0 and 1.', default=0.5)
    parser.add_argument('-o', '--output', help='Stem for the PDB output', default='network')
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
                dist_list.append(1E6)
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
    # count occupancy of each cluster
    clust_occ = [[i, list(clust_ids).count(i)] for i in range(1, n_clusts+1)]
    clust_occ = sorted(clust_occ, key=lambda x:-x[1])
    # renumber the clusters based on this
    old_order = [x[0] for x in clust_occ]
    clust_occs = [x[1] for x in clust_occ]
    # this renumbers the cluster ids for each water so that 1 is the most occupied...
    clust_ids_sorted = np.asarray([old_order.index(x)+1 for x in clust_ids])
    return clust_ids_sorted, clust_occs


def write_average_clusts(filename, clust_wat_ids, n_clusts, n_frames, clust_occs):
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

    Returns
    -------
    centres : dict
        Dictionary containing the average coordinates of each cluster
    """
    # Print cluster occupancy data to screen
    for i in range(1, n_clusts+1):
        print("\tCluster {:2d}:\t{:3d}/{:3d} frames".format(i, len(clust_wat_ids[i]), n_frames))
    print("")
    # Calculate average coordinates and write to PDB
    centres = {}
    with open(filename, 'w') as f:
        f.write("REMARK Average cluster locations written by network_builder.py\n")
        for n in range(1, n_clusts+1):
            # Calculate average coordinates
            av_coords = np.zeros(3)
            for i in clust_wat_ids[n]:
                av_coords += coord_list[i]
            av_coords /= len(clust_wat_ids[n])
            centres[n] = av_coords.copy()
            # Write to file
            f.write('HETATM{0:>5} {1:<4} {2:<4} {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.1f}         {9:>2}  \n'.format(n, 'O00', 'WA1',n, av_coords[0], av_coords[1], av_coords[2], 1., int(clust_occs[n-1]), 'O'))
# calculating disorder of each cluster
    clust_disorder = {}
    for n in range(1, n_clusts+1):
        # Calculate average coordinates
        dists = []
        for i in clust_wat_ids[n]:
            dists.append(np.linalg.norm(coord_list[i] - centres[n]))
        disorder = {}
        disorder['min'] = np.min(np.array(dists))
        disorder['max'] = np.max(np.array(dists))
        disorder['mean'] = np.mean(np.array(dists))
        disorder['stddev'] = np.std(np.array(dists))
        clust_disorder[n] = disorder.copy()
    return centres, clust_disorder


# Function to calcuate a 2d distance matrix between cluster centres

def get_distance_matrix(centres):

    """
    Function to calculate a 2D distance matrix between the centres of water clusters

    Parameters
    ----------
    centres : dict
        Dictionary with keys corresponding to cluster IDs and items corresponding
        to NumPy arrays of cluster centre coordinates

    Returns
    -------
    matrix : np.ndarray
        2D matrix containing the distances between all cluster centres
    """

    matrix = np.zeros((len(centres), len(centres)))
    for i in range(len(centres)):
        for j in range(i+1, len(centres)):
            matrix[i,j] = matrix[j,i] = np.linalg.norm(centres[i+1] - centres[j+1])
    return matrix

def calc_phi(a, b, frame_clust_ids):
    """
    Calculate a correlation score for two water sites as the correlation
    between two binary variables (on/off for each site)

    Parameters
    ----------
    a, b : int
        Cluster IDs of the two water sites
    frame_clust_ids : list
        List of lists containing the cluster IDs present in each frame

    Returns
    -------
    phi : float
        Correlation score (phi) between the two sites
    """
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
            n00 += 1
    phi = (n11*n00 - n01*n10) / ((n11+n01)*(n11+n10)*(n00+n01)*(n00+n10))**0.5
    return phi

def calc_phi_norm(a, b, frame_clust_ids): 
    """ 
    Calculate a correlation score for two water sites as the correlation 
    between two binary variables (on/off for each site) 
 
    Parameters 
    ---------- 
    a, b : int 
        Cluster IDs of the two water sites 
    frame_clust_ids : list 
        List of lists containing the cluster IDs present in each frame 
 
    Returns 
    ------- 
    phi : float 
        Correlation score (phi) between the two sites 
    """ 
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
            n00 += 1
    print n00, n01, n10, n11
    n11 = float(n11)/400.
    n01 = float(n01)/400. 
    n10 = float(n10)/400. 
    n00 = float(n00)/400. 
    na0 = (n00+n01) 
    nb0 = (n00+n10) 
    na1 = (n10+n11) 
    nb1 = (n01+n11) 
    phi = (n11*n00 - n01*n10) / (na0*na1*nb0*nb1)**0.5
    # now need to normalise phi
    if phi >= 0: # need to use phi_max
	if na1 >= nb1:
	    phi_max = (nb1*na0)**0.5/(nb0*na1)**0.5
	    return phi, phi/phi_max
	else:
	    phi_max = (nb0*na1)**0.5/(nb1*na0)**0.5
	    return phi, phi/phi_max
    else: # need to use phi_min
        if na1+nb1 < 1.:
	    phi_min = (nb1*na1)**0.5/(nb0*na0)**0.5
	    return phi, phi/phi_min
	else:
	    phi_min = (nb0*na0)**0.5/(nb1*na1)**0.5
	    return phi, phi/phi_min

def write_pymol(filename, pos, neg):
    """
    Write out a PyMOL script to visualise the system, with lines drawn between correlated waters

    Parameters
    ----------
    filename : str
        Name of output PyMOL script
    pos : list
        List of list describing all positive correlations between clusters
    neg : list
        List of list describing all negative correlations between clusters
    """
    with open(filename, 'w') as f:
        # Write a comment indicating that the code was written by this script
        f.write("# PyMOL visualisation script written by network_builder.py from ProtoMS (www.protoms.org)\n\n")
        # Module loading and initial setup
        f.write("import __main__\n__main__.pymol_argv = [ 'pymol']\n")
        f.write("import pymol\npymol.finish_launching()\nfrom pymol import cmd\nimport os\n\n\n")
        # Read in PDB files
        f.write("# Read in cluster centres\ncmd.load('clusts.pdb')\n")
        f.write("# Read in protein - change this to look for protein_scoop.pdb before using protein.pdb\n")
        f.write("cmd.load('../protein.pdb')\ninputfiles = os.listdir('..')\n\n")
        f.write("for x in inputfiles:\n    if len(x) == 7:\n        if x[3:] == '.pdb':\n")
        f.write("            # This is probably the ligand file\n")
        f.write("            cmd.load('../{}'.format(x))\n            ligname = x[0:3]\n\n")
        # Format visualisation of protein, ligand & waters
        f.write("cmd.hide('lines')  # Remove lines\n\n# Format ligand\ncmd.show('sticks', ligname)\n\n")
        f.write("cmd.hide('sticks','h. and (e. c extend 1)')\n")
        f.write("# Format protein\ncmd.show('cartoon','protein')\ncmd.set('cartoon_transparency',0.5)\n\n")
        f.write("# Format water clusters\ncmd.show('spheres','clusts')\ncmd.set('sphere_scale',0.2)\n")
        f.write("cmd.spectrum('b','blue_red','clusts')\ncmd.label('clusts','resi')\ncmd.orient('clusts')\n\n")
        # Draw correlation lines
        
        f.write("# Draw positive and negative correlations\n")
	if len(pos) > 0:
        	f.write("pos = {}\n".format(pos, neg))
        	f.write("for pair in pos:\n")
        	f.write("    cmd.distance('pos', 'clusts and i. {}'.format(pair[0]), 'clusts and i. {}'.format(pair[1]))\n")
        	f.write("    cmd.color('green','pos')\n\n")
        	f.write("# Hide distance labels on correlation 'bonds'\ncmd.hide('labels','pos')\n")
	if len(neg) > 0:
        	f.write("neg = {}\n\n".format(neg))
        	f.write("for pair in neg:\n")
        	f.write("    cmd.distance('neg', 'clusts and i. {}'.format(pair[0]), 'clusts and i. {}'.format(pair[1]))\n")
        	f.write("    cmd.color('red','neg')\n\n")
        	f.write("# Hide distance labels on correlation 'bonds'\ncmd.hide('labels','neg')\n")
    return None


def get_rep_frames(rep_networks, frame_clust_ids, frame_wat_ids, clust_wat_ids, clust_centres):
    """
    Get representative frames for each of the networks generated

    Parameters
    ----------
    rep_networks : list
        List of lists containing the cluster IDs present in each network
    frame_clust_ids : list
        List of lists containing the cluster IDs present in each frame
    frame_wat_ids : list
        List of lists containing the water observations present in each frame
    clust_wat_ids : dict
        Dictionary containing the water observations contained in each cluster
    clust_centres : dict
        Dictionary containing average positions of each cluster

    Returns
    -------
    network_rep_frames : list
        List of frame representing each of the networks
    """
    # Need to identify which frames contain each network
    network_frame_ids = [[] for net in rep_networks]
    for i, network in enumerate(rep_networks):
        for j, frame in enumerate(frame_clust_ids):
            # Check that the network and frame contain the same clusters - ideally want the whole network
            if len(network) == len(frame) and np.all([wat in frame for wat in network]):
                network_frame_ids[i].append(j)
        # The whole network may not be represented by any frames, so we could look for frames which contain the network
        if len(network_frame_ids[i]) == 0:
            print("\tCan't find an exact frame match for network {}, using frames containing the network".format(i+1))
            for j, frame in enumerate(frame_clust_ids):
                # Checking less strictly now - the rep. network and frame need not be exactly identical
                if np.all([wat in frame for wat in network]):
                    network_frame_ids[i].append(j)
        # If there are still no frames then something is wrong and we raise an exception
        if len(network_frame_ids[i]) == 0:
            raise RuntimeError("Network {} is not present or included in any frame".format(i+1))
    # Now want to find the representative frames
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
                    else:
                        continue
            rmsd = np.sqrt(sum(sq_dists))
            if rmsd < min_rmsd:
                min_rmsd = rmsd
                best_frame = j
        network_rep_frames.append(best_frame)
    return network_rep_frames


if __name__ == "__main__":
    start = time.time()
    # Read command line arguments
    args = get_args()
    
    # Read PDB input data - only water molecules
    reading_time = time.time()  # Measure time spent reading data
    pdbfiles = simulationobjects.PDBSet()
    pdbfiles.read(args.input, resname=args.molecule, skip=args.skip, readmax=9999)
    n_frames = len(pdbfiles.pdbs)
    print('\n{} frames of PDB data read in.\n'.format(n_frames))
    
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
    reading_time = time.time() - reading_time

    # Calculate a 1D distance matrix between all water observations
    dist_time = time.time()
    dist_list = get_distances(frame_wat_ids, coord_list)
    dist_time = time.time() - dist_time
    
    # Cluster waters into discrete hydration sites
    cluster_time = time.time()
    print("Clustering water sites...")
    tree = hierarchy.linkage(dist_list, method=args.linkage)
    clust_ids = hierarchy.fcluster(tree, t=args.distance, criterion='distance')
    n_clusts = max(clust_ids)
    print("\n{} clusters identified:".format(n_clusts))
    cluster_time = time.time() - cluster_time

    # Renumber clusters based on occupancy
    clust_ids_sorted, clust_occs = sort_clusters(clust_ids)

    # Identify which water observations are present in each cluster
    clust_wat_ids = {}  # Store water IDs for each cluster
    for i in range(1, n_clusts+1):
        clust_wat_ids[i] = []
    for i in range(total_wats):
        clust = clust_ids_sorted[i]
        clust_wat_ids[clust].append(i)

    # Prints out the average of each cluster to PDB
    clust_centres,clust_disorder = write_average_clusts("clusts.pdb".format(args.output), clust_wat_ids,
                                         n_clusts, n_frames, clust_occs)
    # Check which clusters are observed in each frame
    clust_frame_ids = {}  # Store frame IDs for each cluster
    frame_clust_ids = [[] for i in range(n_frames)]  # Stores clusters present in each frame
    for i in range(1, n_clusts+1):
        clust_frame_ids[i] = []
        for wat_id in clust_wat_ids[i]:
            for j, frame in enumerate(frame_wat_ids):
                if wat_id in frame:
                    clust_frame_ids[i].append(j)
                    frame_clust_ids[j].append(i)
    # Carry out the correlation analysis
    # We consider waters that are on less than 100% of the time, and more than the cutoff occupancy specified
    correl_time = time.time()
    positives = []  # Store all positives
    pos_close = []  # Store close positives - within 5 A currently
    negatives = []  # Store all negatives
    neg_close = []  # Store close negatives - within 5 A currently
    correlation_water_set = [i for i in range(1, n_clusts+1) if 100 <= clust_occs[i-1] < n_frames] # only considering waters with greater than 25% occupancy
    phidict = {}
    # Starting with pairwise correlations - may figure out n-body later...
    print("Analysing water-water correlations...")
    for x, y in itertools.combinations(correlation_water_set, 2):
        phi_coef = calc_phi(x, y, frame_clust_ids)
        if phi_coef > args.phi:
            print("\tPositive correlation between clusters {:2d} and {:2d} (phi = {:+5.2f})".format(x, y, np.round(phi_coef, 2)))
            positives.append([x, y])
            if np.linalg.norm(clust_centres[x] - clust_centres[y]) < 3.4:
                pos_close.append([x, y])
		phidict[(x,y)] = phi_coef
        elif phi_coef < -args.phi:
            print("\tNegative correlation between clusters {:2d} and {:2d} (phi = {:+5.2f})".format(x, y, np.round(phi_coef, 2)))
            negatives.append([x, y])
            if np.linalg.norm(clust_centres[x] - clust_centres[y]) < 3.4:
                neg_close.append([x, y])
		phidict[(x,y)] = phi_coef
    print("{} positive and {} negative correlations found".format(len(positives), len(negatives)))
    correl_time = time.time() - correl_time
    # Write PyMOL visualisation script - showing close correlations
    write_pymol("pymol-{}.py".format(args.output), pos_close, neg_close)

    distances = get_distance_matrix(clust_centres) 

    all_subnets = []
    all_occupied = [i for i in range(1, n_clusts+1) if 100 <= clust_occs[i-1] ]# this is different to the correlation_water_set as includes 100% occupied water 
    for i in all_occupied: 
	subnet = [i] 
    	for repeat in range(len(all_occupied)**2): # need to loop through a few times
	        for j in all_occupied:
			if j in subnet:
				continue	
			maxdistance = 3.4 + clust_disorder[i]['stddev'] + clust_disorder[j]['stddev'] 
			mindistance = 2.4 - clust_disorder[i]['stddev'] - clust_disorder[j]['stddev']
			if any (mindistance < d < maxdistance for d in [distances[x-1][j-1] for x in subnet]):
				if all (d > mindistance for d in [distances[x-1][j-1] for x in subnet]): # this prevents short range clashes
					pair = tuple(sorted([i,j]))
					if pair in phidict.keys():
						if phidict[pair] > -args.phi: #this prevents negative correlations
							subnet.append(j)
					else:
						subnet.append(j)
	if len(subnet) > 1: # dont want single water subnets
		all_subnets.append(sorted(subnet))
    print 'SUBNETWORKS'
    subnetworks = set(tuple(row) for row in all_subnets)
    subnetworks = [list(row) for row in subnetworks]
    subnetworks.sort(key=len)
    for a in subnetworks:
	print a
    print 'THERE ARE ',len(subnetworks),'SUBNETWORKS'
    print 'now for filtering down the number of subnetworks'

    quit()
    print("{} network(s) built.".format(len(networks)))

    # Generate a representative frame for each network
    print("\nFinding representative frames...")
    network_rep_frames = get_rep_frames(networks, frame_clust_ids, frame_wat_ids,
                                        clust_wat_ids, clust_centres)
    network_time = time.time() - network_time
   
    # Write out representative frames to PDB files..
    writing_time = time.time()
    print("\nWriting PDB output...")
    for i, j in enumerate(network_rep_frames):
        # Write out water only for rep. frames
        pdbfiles.pdbs[j].write("{}-{:02d}-wat.pdb".format(args.output, i+1))
        # Read in whole system - only want to read a single frame to save time & memory
        pdbfiles2 = simulationobjects.PDBSet()
        pdbfiles2.read(args.input, skip=args.skip+j, readmax=1)
        pdbfiles2.pdbs[0].write("{}-{:02d}-all.pdb".format(args.output, i+1))
    writing_time = time.time() - writing_time

    # Print out timings of analysis - for testing and reference
    print("\nNetwork analysis completed in {} seconds".format(np.round(time.time()-start, 1)))
    print("\tReading data:          {:7.1f} seconds".format(np.round(reading_time, 1)))
    print("\tDistance calculations: {:7.1f} seconds".format(np.round(dist_time, 1)))
    print("\tClustering:            {:7.1f} seconds".format(np.round(cluster_time, 1)))
    print("\tCorrelation analysis:  {:7.1f} seconds".format(np.round(correl_time, 1)))
    print("\tNetwork construction:  {:7.1f} seconds".format(np.round(network_time, 1)))
    print("\tWriting data:          {:7.1f} seconds\n".format(np.round(writing_time, 1)))