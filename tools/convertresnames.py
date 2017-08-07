"""
Script to convert residue names to ProtoMS format, before applying the atom name conversion
Some residues, such as HIE must be formally labelled as such or they will not be read correctly.
This script seeks to automate this feature of the setup.

Can be executed alone or used to import functions

Marley Samways
"""

import logging

logger = logging.getLogger('protoms')


def rename_residues(pdb_in):
    """
    Function to read in all residues, identify those which need changing, and rename them appropriately
    """
    pdb_out = pdb_in.copy()
    # Count histidine conversions
    n_his_in = 0
    n_hie_out = 0
    n_hip_out = 0
    # Count cysteine conversions
    n_cys_in = 0
    n_cyx_out = 0
    # Loop over residues and check for atoms which are indicative of a different form...
    for resnum in pdb_in.residues:
        residue = pdb_in.residues[resnum]
        if residue.name.upper() == "HIS":
            n_his_in += 1
            epsilon = False  # True if NE is protonated
            delta = False    # True if ND is protonated
            for atomnum in range(len(residue.atoms)):
                atomname = residue.atoms[atomnum].name.upper()
                if atomname == "HND" or atomname == "HD1":
                    delta = True
                elif atomname == "HNE" or atomname == "HE2":
                    epsilon = True
            if epsilon:
                if delta:
                    # Double protonated
                    n_hip_out += 1
                    pdb_out.residues[resnum].name = "HIP"
                else:
                    # Protonated only at NE
                    n_hie_out += 1
                    pdb_out.residues[resnum].name = "HIE"
            else:
                # Protonated only at ND => leave name as HIS
                continue
        elif residue.name.upper() == "CYS":
            n_cys_in += 1
            bridging = True  # Assume the residue is bridging, unless an S-H is found
            for atomnum in range(len(residue.atoms)):
                atomname = residue.atoms[atomnum].name.upper()
                if atomname == "HG":
                    bridging = False
            if bridging:
                # Bridging => rename to CYX
                n_cyx_out += 1
                pdb_out.residues[resnum].name = "CYX"
            else:
                # Not bridging => leave name as CYS
                continue
    logger.info("%i/%i HIS residues were renamed to HIP"%(n_hip_out, n_his_in))
    logger.info("%i/%i HIS residues were renamed to HIE"%(n_hie_out, n_his_in))
    logger.info("%i/%i CYS residues were renamed to CYX"%(n_cyx_out, n_cys_in))
    return pdb_out


def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="Script to convert residue names in a protein pdb-file to ProtoMS style")
    parser.add_argument('-p','--protein',help="the protein PDB-file")
    parser.add_argument('-o','--out',help="the output PDB-file",default="protein_rpms.pdb")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    import simulationobjects
    
    args = get_args()
    logger = simulationobjects.setup_logger("convertresnames_py.log")
    
    protein_in = simulationobjects.PDBFile(filename=args.protein)
    protein_out = rename_residues(protein_in)
    protein_out.write(args.out)





