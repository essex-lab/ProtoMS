# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" Routines to convert atom names in a PDB-file

This module defines two public functions:
read_convfile
pdb2pms
    
Can be executed from the command line as a stand-alone program
"""

import logging

logger = logging.getLogger('protoms')

def read_convfile(file=None, inmode=None,outmode=None):
    """ 
    Reads a conversion file into a dictionary

    Parameters
    ----------        
    file : string
    	the filename of the conversion file
    inmode : string 
    	the style of the input PDB-file
    outmode : string 
    	the style of the output PDB-file

    Returns
    -------
    a dictionary of the conversion
    
    Raises
    ------
    ValueError
    	if any of the parameters are none
    """

    if file is None:
        raise ValueError("You must pass file!")
    if inmode is None:
        raise ValueError("You must specify a mode!")
    if outmode is None:
        raise ValueError("You must specify a mode!")    
    conv = {}
    stream = open(file,'r')
    buffer = stream.readlines()
    stream.close()
    for line in buffer:
        if line.startswith('backbone'):
            key = 'backbone'
            conv[key] = {}
            continue
        if line.startswith('residue'):
            elems = line.split()
            key = elems[1]
            conv[key] = {}
            continue        
        if line.startswith('atom'):
            elems = line.split()
            # Locate inmode
            found = False
            for x in range(0,len(elems)):
                if elems[x] == inmode:
                    found = True
                    old = elems[x+1]
                    break
            if not found:
                continue
            #Now locate outmode
            found = False
            for x in range(0,len(elems)):
                if elems[x] == outmode:
                    found = True
                    new = elems[x+1]
                    break
            if not found:
                continue            
                #raise "Mode %s not present for %s " % (mode,line)
            conv[key][old] = new
    return conv

def _checkpdb(pdb_in,conversion) :
    """
    Checks if all atom names are already in the PDB file

    Parameters
    ----------        
    pdb_in : PDBFile 
      	the pdb file
    conversion :
    	a dictionary with conversion

    Returns
    -------
    bool 
      True if all atom names are already in the PDB file
      False if not
    """

    for resnum in pdb_in.residues:
        residue = pdb_in.residues[resnum]

        if residue.name.upper() == "HID": # Valid only when the pdb input file is from AMBER.
            resname = "HIS"
            pdb_out.residues[resnum].name = resname
        else:
            resname = residue.name.upper()

        for atomnum in range(len(residue.atoms)):
            atomname = residue.atoms[atomnum].name.upper()
            # Check if the atom name is either in the backbone or a sidechain residue
            if not (atomname in conversion["backbone"].values() or atomname in conversion[resname].values()) :
              return False    
    return True

def pdb2pms(pdb_in,forcefield,conversion_file):
    """ 
    Convert atom names to ProtoMS standard

    The protein object is not modified by this routine, but the Residue and Atom objects are.

    Parameters
    ----------        
    pdb_in : PDBFile 
      	the pdb file to modify inline
    forcefield : string 
    	the style of the input PDB-file
    conversion_file :
    	a file with conversion rules

    Returns
    -------
    a PDBFile instance with ProtoMS names
    """

    logger.debug("Running pdb2pms with arguments: ")
    logger.debug("\tpdb_in          = %s"%pdb_in) 
    logger.debug("\tforcefield      = %s"%forcefield) 
    logger.debug("\tconversion_file = %s"%conversion_file) 
    logger.debug("This will rename atoms in a PDB-file to match ProtoMS naming convention")

    pdb_out = pdb_in.copy()
    conversion = read_convfile(conversion_file,inmode=forcefield,outmode="protoms")
    if _checkpdb(pdb_in,conversion) :
      logger.info("It seems that %s already contains the correct ProtoMS names. Aborting the conversion."%pdb_in.name)
      return pdb_in  
    for resnum in pdb_in.residues:
        residue = pdb_in.residues[resnum]
        if residue.name.upper() == "HID":									# Valid only when the pdb input file is from AMBER.
            resname = "HIS"
            pdb_out.residues[resnum].name = resname
        else:
            resname = residue.name.upper()
        for atomnum in range(len(residue.atoms)):
            atomname = residue.atoms[atomnum].name.upper()
            if atomname in conversion["backbone"]:
                newName = conversion["backbone"][atomname]
            elif atomname in conversion[resname]:
                newName =  conversion[resname][atomname]
            elif atomname[1:len(atomname)] + atomname[0] in conversion[resname]:
                newName = conversion[resname][(atomname[1:len(atomname)] + atomname[0])]
            elif atomname[-1] + atomname[0:(len(atomname)-1)] in conversion[resname]:
                newName = conversion[resname][(atomname[-1] + atomname[0:(len(atomname)-1)])]
            else:
                newName = atomname
                logger.warning("Warning: atom %s in residue %s %i not found in %s. This atom has not been converted. Check atom names are valid."  %(atomname, resname, resnum, conversion_file))
            pdb_out.residues[resnum].atoms[atomnum].name = newName
            pdb_out.residues[resnum].atoms[atomnum].resname = resname 
    return pdb_out

def get_arg_parser():
    
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program convert atom names in a protein pdb-file to ProtoMS style")
    parser.add_argument('-p','--protein',help="the protein PDB-file")
    parser.add_argument('-o','--out',help="the output PDB-file",default="protein_pms.pdb")
    parser.add_argument('-s','--style',help="the style of the input PDB-file",default="amber")
    parser.add_argument('-c','--conversionfile',help="the name of the file with conversion rules",default="atomnamesmap.dat")
    return parser

if __name__ == "__main__":

    import simulationobjects

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("convertatomnames_py.log")

    protein = simulationobjects.PDBFile(filename=args.protein)
    protein_out = pdb2pms(protein,args.style,args.conversionfile)
    protein_out.write(args.out)
