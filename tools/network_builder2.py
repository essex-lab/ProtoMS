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


def first_solvation_shell(ligand, clust_centers, clust_disorders):

	lig_heavy = []

	with open(ligand, 'r') as f:
		for line in f:
			if line.split()[0] == 'HETATM':
				if line.split()[2][0] != 'H':
				    lig_heavy.append([float(line[27:38]),float(line[39:47]),float(line[47:55])])
	lig_heavy = np.asarray(lig_heavy)
	first_shell_results = []

	for i in range(1, len(clust_centers)+1):
		first_shell = False
		for atom in lig_heavy:
			distance = np.linalg.norm(clust_centers[i] - atom)
			cutoff = 3.4+clust_disorders[i-1]
			if distance <= cutoff :
				first_shell = True
		first_shell_results. append(first_shell)
	in_first_shell = [x for x,res in enumerate(first_shell_results,1) if res == True]
	return in_first_shell


def get_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser('network_builder.py')
    parser.add_argument('-i', '--input', default='all.pdb',
                        help='PDB file containing input frames. Currently accepts only one file')
    parser.add_argument('-m', '--molecule', default='WA1',
                        help='Residue name of water molecules')
    parser.add_argument('-a', '--atom', default='O00',
                        help='Name of atom to take as molecule coordinates')
    parser.add_argument('-s', '--skip', type=int, default=0,
                        help='Number of frames to skip')
    parser.add_argument('-d', '--distance', type=float, default=3.0,
                        help='Distance cutoff for clustering. Default=3.0 Angs')
    parser.add_argument('--linkage', default='average',
                        help='Linkage method for hierarchical clustering')
    parser.add_argument('-o', '--occupancy', type=float, default=0.0,
                        help='Occupancy cutoffs for waters to be considered. Between 0 and 1.')
    parser.add_argument('-c', '--correlation', type=float, default=0.1,
                        help='Percentage correlation (greater than random) deemed significant. Between 0 and 1.')
    parser.add_argument('-l', '--ligand',
                        help='Input PDB file of the ligand. Needed for the first solvation shell')
    parser.add_argument('--firstshell', action='store_true', default=False,
                        help='Only perform network analysis on first solvation shell')
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
    dist_mat : np.ndarray
        2D distance matrix
    """
    print("Calculating water-water distances...")
    n_wats = len(coord_list)
    dist_list = []
    dist_mat = np.zeros((n_wats, n_wats))
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
                dist = np.linalg.norm(coord_list[i]-coord_list[j])
                dist_mat[i,j] = dist_mat[j,i] = dist
            else:
                dist = np.linalg.norm(coord_list[i]-coord_list[j])
                dist_list.append(dist)
                dist_mat[i,j] = dist_mat[j,i] = dist
                frame_last = False  # Stops checking if subsequent j values are in the same frame
    return dist_list, dist_mat


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

def write_average_clusts(filename, first_shell, centers, clust_occs):
    with open(filename, 'w') as f:
        f.write("REMARK Average cluster locations written by network_builder.py\n")
        for n in first_shell:
            # Write to file
            f.write('HETATM{0:>5} {1:<4} {2:<4} {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.1f}         {9:>2}  \n'.format(n, 'O00', 'WA1',n, centers[n][0], centers[n][1], centers[n][2], 1., int(clust_occs[n-1]), 'O'))

    return

def calc_average_clusts(clust_wat_ids, n_clusts, n_frames, clust_occs):
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
    for i in range(1, n_clusts +1):
        print("\tCluster {:2d}:\t{:3d}/{:3d} frames".format(i, len(clust_wat_ids[i]), n_frames))
    print("")
    # Calculate average coordinates and write to PDB
    centres = {}
    for n in range(1, n_clusts+1):
        # Calculate average coordinates
        av_coords = np.zeros(3)
        for i in clust_wat_ids[n]:
            av_coords += coord_list[i]
        av_coords /= len(clust_wat_ids[n])
        # Find real position closest to cluster centre
        min_dist = 1E6
        for i in clust_wat_ids[n]:
            dist = np.linalg.norm(coord_list[i]-av_coords)
            if dist < min_dist:
                centres[n] = coord_list[i].copy()
                min_dist = dist
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


def get_cluster_distances(frame_clust_ids, frame_wat_ids, clust_wat_ids, orig_matrix, centres):

    """
    Function to calculate a 2D distance matrix between the observations contained in each cluster

    Parameters
    ----------
    frame_clust_ids : list
        List of which clusters are present in each frame
    frame_wat_ids : list
        List of which water IDs are present in each frame
    clust_wat_ids : dict
        List of the water IDs contained by each cluster
    coord_list : list
        List of water oxygen coordinates

    Returns
    -------
    matrix : np.ndarray
        2D matrix containing the distances between all cluster centres
    """
    n_clusts = len(clust_frame_ids)
    matrix = np.zeros((n_clusts, n_clusts))
    for i in range(n_clusts):
        for j in range(i+1, n_clusts):
            dists = []  # List of observed distances
            for wat_i in clust_wat_ids[i+1]:
                for frame in frame_wat_ids:
                    if wat_i not in frame:
                        continue
                    for wat_j in clust_wat_ids[j+1]:
                        if wat_j not in frame:
                            continue
                        dists.append(orig_matrix[wat_i, wat_j])
                        break
            # If clusters are never observed together, find the distance between the means 
            if len(dists) == 0:
                matrix[i,j] = matrix[j,i] =  -1
            else:
                matrix[i,j] = matrix[j,i] = np.mean(dists)
    return matrix


def calc_correlation(a, b, frame_clust_ids):
    a_on = 0.0
    b_on = 0.0
    both = 0.0
    n_frames = 0.0
    for frame in frame_clust_ids:
        n_frames += 1
        if a in frame:
            a_on += 1
        if b in frame:
            b_on += 1
        if a in frame and b in frame:
            both += 1
    a_on /= n_frames
    b_on /= n_frames
    both /= n_frames
    random = a_on * b_on 
    correlation = both - random
    return correlation


def tanimoto(list_a, list_b):
    a = float(len(list_a))
    b = float(len(list_b))
    c = float(len([x for x in list_a if x in list_b]))
    S = c / (a + b - c)
    return S


def write_pymol(filename, pos, neg, nets):
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

	f.write("# Draw networks\n")
	colors = {1:'yellow', 2:'blue', 3:'orange',4:'purple',5:'cyan',6:'green',7:'green',8:'green',9:'green',10:'green'}
	for i, net in enumerate(nets,1):
		f.write("net"+str(i)+" = {}\n\n".format(net))
		f.write("for pair in net"+str(i)+":\n")
		f.write("    cmd.distance('net"+str(i)+"', 'clusts and i. {}'.format(pair[0]), 'clusts and i. {}'.format(pair[1]))\n")
		#f.write("    cmd.color('"+colors[i]+"','net"+str(i)+"')\n\n")
 	f.write("# Hide distance labels on network 'bonds'\ncmd.hide('labels','net*')\n")

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
    # Read command line arguments
    start = time.time()
    args = get_args()

    # Print out what type of analysis will be performed
    if args.firstshell:
        if args.ligand is None:
            raise RuntimeError("\nCan't carry out first solvation shell analysis without a ligand file!\n")
        else:
            print("\nOnly analysing waters in the first solvation shell...")
    else:
        print("\nCarrying out network analysis on all GCMC water sites...")
        if args.ligand is not None:
            print("Ligand file {} will be ignored".format(args.ligand))
            args.ligand = None
    
    # 
    # Read in data
    # 

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

    # 
    # Calculate water-water distances and cluster hydration sites
    # 

    # Calculate a 1D distance matrix between all water observations
    dist_time = time.time()
    dist_list, dist_matrix = get_distances(frame_wat_ids, coord_list)
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
    clust_centres, clust_disorder = calc_average_clusts(clust_wat_ids, n_clusts, n_frames, clust_occs)

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

    # Calculating cluster-cluster distances
    clust_dist_time = time.time()
    distances = get_cluster_distances(frame_clust_ids, frame_wat_ids, clust_wat_ids, dist_matrix, clust_centres) # mean framewise distance between waters in cluster
    clust_dist_time = time.time() - clust_dist_time

    # 
    # Analyse correlations between clusters
    # 

    # Carry out the correlation analysis
    # We consider waters that are on less than 100% of the time, and more than the cutoff occupancy specified
    correl_time = time.time()
    positives = []  # Store all positives
    negatives = []  # Store all negatives
    correlation_water_set = [i for i in range(1, n_clusts+1) if args.occupancy*n_frames <= clust_occs[i-1] < n_frames]  # only considering waters with greater than specified occupancy
    correl_dict = {}
    # Starting with pairwise correlations - may figure out n-body later...
    print("Analysing water-water correlations...")
    for x, y in itertools.combinations(correlation_water_set, 2):
	if 2.4 <= distances[x-1][y-1] <= 3.4:
            correl = calc_correlation(x, y, frame_clust_ids)
            if correl > args.correlation:
                print("\tPositive correlation between clusters {:2d} and {:2d} ({:+6.2f})".format(x, y, np.round(100*correl, 2)))
                positives.append([x, y])
    		correl_dict[(x,y)] = correl
            elif correl < -args.correlation:
                print("\tNegative correlation between clusters {:2d} and {:2d} ({:+6.2f})".format(x, y, np.round(100*correl, 2)))
                negatives.append([x, y])
    		correl_dict[(x,y)] = correl
    print("{} positive and {} negative correlations found".format(len(positives), len(negatives)))
    correl_time = time.time() - correl_time

    # 
    # Build networks based on hydrogen bonds
    # 

    # Start building networks    
    network_time = time.time()
    # Separate sites into first shell, if requested
    if args.ligand is not None:	
        net_sites = first_solvation_shell(args.ligand, clust_centres, disorder_list)
        if len(net_sites) == 0:
            raise RuntimeError("No waters are in the first solvation shell")
        print("Clusters in the first solvation shell:")
        print("\t{}".format(net_sites))
    else:
        net_sites = [i for i in range(1, n_clusts+1)]

    # Construct networks
    networks = []
    for i in net_sites:
        # Ignore sites with too low occupancy
        if args.occupancy * n_frames > clust_occs[i-1]:
            continue
        sub_net = [i]  # Start with one water
        # Iteratively add waters to the network
        while True:
            building = False  # Mark whether a site has been added to the network
            for j in net_sites:
                # Ignore if j has already been counted or is too low in occupancy
                if j in sub_net or args.occupancy * n_frames > clust_occs[j-1]:
                    continue
                # Make sure that water j is H-bonding to at least one water in the network
                if not any(2.4 <= d <= 3.4 for d in [distances[x-1, j-1] for x in sub_net]):
                    continue
                # Make sure that j does not clash with any waters in the network
                if any(d < 2.4 for d in [distances[x-1, j-1] for x in sub_net]):
                    continue
                # Check that there is at least one frame where all of these waters have been observed together
                if not any((all(x in frame for x in sub_net) and j in frame) for frame in frame_clust_ids):
                    continue
                # Check if there are any negative correlations between j and anything in the network
                clashing = False
                for x in sub_net:
                    pair = tuple(sorted([x, j]))
                    if pair in correl_dict.keys():
                        if correl_dict[pair] < -args.correlation:
                            clashing = True
                # If there is no clash, the water can join the network
                if not clashing:
                    sub_net.append(j)
                    building = True
            # Break the loop if no waters have been added
            if not building:
                break
        # Add to the list, if a network (N > 1) has been formed
        if len(sub_net) > 1:
            networks.append(sorted(sub_net))
    # Remove any repeated networks and sort by size
    networks = set(tuple(net) for net in networks)
    networks = [list(net) for net in networks]
    networks.sort(key=len)
    print("\n{} network(s) built.".format(len(networks)))

    # Calculate occupancies of these networks
    net_occs = []
    for i in range(len(networks)):
        network_occ = 0
        for frame in frame_clust_ids:
            if all(wat in frame for wat in networks[i]):
                network_occ += 1
        net_occs.append(network_occ)

    # Filter networks based on the extent of their overlap
    remove_networks = []  # Networks to remove
    network_comparisons = []  # Record comparisons made
    while True:
        # First want to find the most similar pair of networks which have not been compared
        max_similarity = -1
        for i in range(len(networks)):
            if i in remove_networks:
                continue
            for j in range(i+1, len(networks)):
                if j in remove_networks:
                    continue
                if [i, j] in network_comparisons or [j, i] in network_comparisons:
                    continue
                similarity = tanimoto(networks[i], networks[j])
                if similarity > max_similarity:
                    max_similarity = similarity
                    max_i, max_j = i, j
        # Check if the highest similarity is greater than some threshold
        if max_similarity < 0.5:
            break
        # First prioritise most common network
        if net_occs[max_i] > net_occs[max_j]:
            remove_networks.append(max_j)
        elif net_occs[max_j] > net_occs[max_i]:
            remove_networks.append(max_i)
        # Then prioritise the larger network
        elif len(networks[max_i]) > len(networks[max_j]):
            remove_networks.append(max_j)
        elif len(networks[max_j]) > len(networks[max_i]):
            remove_networks.append(max_i)
        # If they cannot be distinguished, then remove neither & note the comparison
        network_comparisons.append([max_i, max_j])
    filtered_networks = [networks[i] for i in range(len(networks)) if i not in remove_networks]
    print("Filtered to {} networks".format(len(filtered_networks)))

    # Print out similarities of networks
    print("\nNetwork similarities:")
    for i in range(len(filtered_networks)):
        for j in range(i+1, len(filtered_networks)):
            S = tanimoto(filtered_networks[i], filtered_networks[j])
            print("\t{:2d} & {:2d} : S = {:.3f}".format(i+1, j+1, np.round(S, 3)))

    # Generate a representative frame for each network
    print("\nFinding representative frames...")
    network_rep_frames = get_rep_frames(filtered_networks, frame_clust_ids, frame_wat_ids,
                                        clust_wat_ids, clust_centres)
    network_time = time.time() - network_time

    # 
    # Write network and cluster data to PDB files and a PyMOL script
    # 
   
    # Determine H-bonding pairs in networks (for visualisation)
    all_pairs = []
    for net in filtered_networks:
	pairs = []
	for i,j in itertools.combinations(net,2):
	    if 2.4 <= distances[i-1][j-1] <= 3.4:
                pairs.append([i,j])
	all_pairs.append(pairs)

    # Write clusters to PDB file and also a PyMOL script for visualisation
    write_average_clusts("clusts.pdb".format(args.output), net_sites, clust_centres, clust_occs)
    write_pymol("pymol-{}.py".format(args.output), positives, negatives,all_pairs)

    # Write out representative frames to PDB files..
    writing_time = time.time()
    print("\nWriting PDB output...")
    for net_id, frame in enumerate(network_rep_frames):
        # Get the representative waters for this frame
        net_wats = []
        for clust_id in filtered_networks[net_id]:
            for wat_id in clust_wat_ids[clust_id]:
                if wat_id in frame_wat_ids[frame]:
                    net_wats.append(wat_id)
        # Write out waters only for rep. frame
        with open("{}-{:02d}-wat.pdb".format(args.output, net_id+1), 'w') as f:
            for n, wat_id in zip(filtered_networks[net_id], net_wats):
                for i, atom in enumerate(wat_list[wat_id].atoms):
                    atom_id = (n - 1) * len(wat_list[wat_id].atoms) + i + 1
                    f.write('ATOM  {0:>5} {1:<4} {2:<4} {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.2f}\n'.format(atom_id, atom.name, 'WA1', n, atom.coords[0], atom.coords[1], atom.coords[2], clust_occs[n-1]/float(n_frames), clust_occs[n-1]))
                f.write('TER\n')
            f.write('END')
        # Read in whole system - only want to read a single frame to save time & memory
        pdbfiles2 = simulationobjects.PDBSet()
        pdbfiles2.read(args.input, skip=args.skip+frame, readmax=1)
        pdbfiles2.pdbs[0].write("{}-{:02d}-all.pdb".format(args.output, net_id+1))
    writing_time = time.time() - writing_time

    # Print out timings of analysis - for testing and reference
    print("\nNetwork analysis completed in {} seconds".format(np.round(time.time()-start, 1)))
    print("\tReading data:          {:7.1f} seconds".format(np.round(reading_time, 1)))
    print("\tDistance calculations: {:7.1f} seconds".format(np.round(dist_time, 1)))
    print("\tClustering:            {:7.1f} seconds".format(np.round(cluster_time, 1)))
    print("\tCluster distances:     {:7.1f} seconds".format(np.round(clust_dist_time, 1)))
    print("\tCorrelation analysis:  {:7.1f} seconds".format(np.round(correl_time, 1)))
    print("\tNetwork construction:  {:7.1f} seconds".format(np.round(network_time, 1)))
    print("\tWriting data:          {:7.1f} seconds\n".format(np.round(writing_time, 1)))

