# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

""" 
Routines to build a protein scoop

This module defines a single public function
scoop
    
Can be executed from the command line as a stand-alone program
"""

import os
import operator
import itertools

import numpy as np

import simulationobjects

def scoop ( protein, ligand, out_file = '',
            innercut = 16, outercut  = 20, 
            flexin = 'full', flexout = 'sidechain', 
            excluded = [], added = [] ):

    """Function to generate scoop from protein structure

    scoop ( protein,  ligand, out_file = '', 
            innercut = 16, outercut = 20,
            flexin = 'full', flexout = 'sidechain',
            excluded = [], added = [] )

    Parameters
    ----------
    protein : simulationobjects.pdb
        pdb instance for the protein to be scooped
    ligand : simulationobjects.pdb or string
        Either pdb instance for ligand to be scooped around, or string giving
        the name of a file containing 3 floats to act as coords for scoop centre,
        or a string with three floating point numbers
    out_file : string, optional
        File name to write scooped protein. Default of '' does not write a file
    innercut : float, optional
        Maximum distance from ligand defining inner region of the scoop
    outercut : float, optional
        Maximum distance from ligand defining outer region of the scoop
    flexin : string, optional
        Gives the degree of flexibility for residues of the inner region
        Can be 'rigid', 'sidechain' or 'flexible'
    flexout : string, optiniol
        As flexin but for residues of the outer scoop
    excluded : list of integers
        List of indices for residues to be excluded from scoop
    added : list of integers
        List of indices for residues to be included in outer scoop

    Returns
    -------
    A simulationobjects.PDBFile instance representing the scooped protein """

    centerarray = None
    if isinstance(ligand,basestring) :
        if os.path.isfile(ligand) :
          centerarray = np.loadtxt ( ligand )
        else :
          centerarray = np.array(ligand.strip().split(),dtype=float)

    # Either if a string or a file with center coordinates was passed, 
    # we will make it into a dummy pdb object
    if centerarray is not None :
        ligand = simulationobjects.PDBFile ()
        ligand.residues = { 0: simulationobjects.Residue () }
        ligand.residues[0].addAtom ( simulationobjects.Atom ( coords = centerarray ) )
       
    assert flexin in [ 'sidechain', 'flexible', 'rigid' ]
    assert flexout in [ 'sidechain', 'flexible', 'rigid' ]

    if len ( ligand.residues ) > 1:
        print "More than one ligand in input. Scooping around everything..."
    
    #Build inner and outer lists
    #If any heavy atom of a residue falls within the cutoff distance of
    #any atom of the ligand, add then to the appropriate list
    outer = []
    inner = []
    for res in protein.residues:
        kill = True
        in_kill = True
        for atom in protein.residues[res].atoms:
            if ( atom.name.startswith( ( 'H', 'h', '1' ) )
                 or atom.name in ( 'N', 'C', 'O' ) ):
                continue
            for lig in ligand.residues.itervalues():
                for lat in lig.atoms:
                    if lat.name.startswith ( ( 'H', 'h' ) ):
                        continue
                    distance = np.linalg.norm ( atom.coords - lat.coords )
                    if distance < outercut:
                        kill = False
                    if distance < innercut:
                        in_kill = False

        if not kill and in_kill and not res in excluded:
            outer.append ( res )
        if kill and res in added:
            outer.append ( res )
        if not in_kill:
            inner.append ( res )

    # Also eliminate Xray waters beyond outer cutoff
    waters = []
    for mol in protein.solvents:
        kill = True
        for atom in protein.solvents[mol].atoms:
            for lig in ligand.residues.itervalues():
                for lat in lig.atoms:
                    distance = np.linalg.norm ( atom.coords - lat.coords )
                    if distance < outercut:
                        kill = False
                        break
        if not kill:
            waters.append(mol)


    both = sorted ( inner + outer )

    #All CYS must be fixed to preserve disulphide bridges
    #this may be improved in future
    rigid = [ res for res in both if protein.residues[res].name == 'CYS' ]
    backBoneRigid = []

    if flexout in [ 'rigid', 'sidechain' ]:
        backBoneRigid += [ res for res in outer ]
    if flexout == 'rigid':
        rigid += [ res for res in outer ]
    # Same thing for inner residues
    if flexin in [ 'rigid', 'sidechain' ]:
        backBoneRigid += [ res for res in inner ]
    if flexin == 'rigid':
        rigid += [ res for res in inner ]

    outres = []
    rigidBB = []
    rigidSC = []
    count = 0
    for res in both:
        outres.append(protein.residues[res])
        count += 1
        if res in backBoneRigid:
            rigidBB.append(count)
        if res in rigid:
            rigidSC.append(count)

    # Need to turn residue lists into nice string of residue ranges for ProtoMS to interpret
    # Not crazy about this but works. 
    outBB = ''
    for k, g in itertools.groupby ( enumerate ( rigidBB ), key = lambda ( i, x ): i - x ):
        r = map ( operator.itemgetter ( 1 ), g )
        if len ( r ) > 1:
            outBB += '%d-%d, ' % ( min ( r ), max ( r ) )
        else:
            outBB += '%d, ' % r[0]
    outBB = outBB[:-2]

    outSC = ''
    for k, g in itertools.groupby ( enumerate ( rigidSC ), key = lambda ( i, x ): i - x ):
        r = map ( operator.itemgetter ( 1 ), g )
        if len ( r ) > 1:
            outSC += '%d-%d, ' % ( min ( r ), max ( r ) )
        else:
            outSC += '%d, ' % r[0]
    outSC = outSC[:-2]

    #Purge residues outside the outer scoop from the protein pdb and save it
    for res in protein.residues.keys():
        if res not in both + waters:
            protein.residues.pop ( res )


    header  = "REMARK Scoop of %s\n" % protein
    header += "REMARK Inner Region : %8.2f Angstrom radius\n" % innercut
    header += "REMARK Outer Region : %8.2f Angstrom radius\n" % outercut
    header += "REMARK Num Residues %d ( %d inner %d outer )\n" % (len(inner)+len(outer),
                                                             len(inner),len(outer))
    header += "REMARK %d residues have a fixed backbone\n" % (len(rigidBB))
    header += "REMARK %d residues are fixed\n" % (len(rigidSC))
    header += "REMARK flexibility of the inner part : %s\n" % flexin
    header += "REMARK flexibility of the outer part : %s\n" % flexout
    header += "REMARK ProtoMS keyword to use\n"
    header += "REMARK chunk fixbackbone 1 %s\n" % outBB
    header += "REMARK chunk fixresidues 1 %s\n" % outSC

    header += "REMARK Xray Water within %8.2f Angstrom\n" % outercut
    header += "REMARK of the ligand\n"
    
    protein.header = header

    if out_file:
        protein.write ( out_file, header = header,renumber=True )

    return protein

if __name__ == "__main__":

    import argparse

    import simulationobjects

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program scoop a protein pdb-file")
    parser.add_argument('-p','--protein',help="the protein PDB-file")
    parser.add_argument('-l','--ligand',help="the ligand PDB-file")
    parser.add_argument('-o','--out',help="the output PDB-file",default="scoop.pdb")
    parser.add_argument('--center',help="the center of the scoop, if ligand is not available, either a string or a file with the coordinates",default="0.0 0.0 0.0")
    parser.add_argument('--innercut',type=float,help="maximum distance from ligand defining inner region of the scoop",default=16.0)
    parser.add_argument('--outercut',type=float,help="maximum distance from ligand defining outer region of the scoop",default=20.0)
    parser.add_argument('--flexin',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the inner region",default="flexible")
    parser.add_argument('--flexout',choices=[ 'sidechain', 'flexible', 'rigid' ],help="the flexibility of the inner region",default="sidechain")
    parser.add_argument('--excluded',nargs="+",type=int,help="a list of indices for residues to be excluded from scoops",default=[])
    parser.add_argument('--added',nargs="+",type=int,help="a list of indices for residues to be included in outer scoops",default=[])
    args = parser.parse_args()

    if args.ligand is None :
      ligand = args.center
    else :
      ligand = simulationobjects.PDBFile(filename=args.ligand)
    protein = simulationobjects.PDBFile(filename=args.protein)
    scoop(protein,ligand,args.out,args.innercut,args.outercut,args.flexin,args.flexout,args.excluded,args.added)
  
