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


def pdb2pms(pdb_in,forcefield,conversion_file):
    """ 
    Convert atom names to ProtoMS standard

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

    pdb_out = pdb_in
    conversion = read_convfile(conversion_file,inmode=forcefield,outmode="protoms")
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
                print "Warning: atom %s in residue %s %i not found in %s. This atom has not been converted. Check atom names are valid."  %(atomname, resname, resnum, conversion_file)
            pdb_out.residues[resnum].atoms[atomnum].name = newName
            pdb_out.residues[resnum].atoms[atomnum].resname = resname 
    return pdb_out


if __name__ == "__main__":

    import argparse
    import simulationobjects

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program convert atom names in a protein pdb-file to ProtoMS style")
    parser.add_argument('-p','--protein',help="the protein PDB-file")
    parser.add_argument('-o','--out',help="the output PDB-file",default="protein_pms.pdb")
    parser.add_argument('-s','--style',help="the style of the input PDB-file",default="amber")
    parser.add_argument('-c','--conversionfile',help="the name of the file with conversion rules",default="atomnamesmap.dat")
    args = parser.parse_args()

    protein = simulationobjects.PDBFile(filename=args.protein)
    protein_out = pdb2pms(protein,args.style,args.conversionfile)
    protein_out.write(args.out)
