#!/usr/bin/env python
# Julien Michel November 2004
# 1) Read a connectivity template file in memory
# 2) Read a parameter template file in memory
# 3) Assign parameters to the templates
# 4) Write a template section of a forcefield
#
import os,sys,getopt,copy

class Residue(object):
    def __init__(self,name='',backbone=None,rotate=0.0,translate=0.0,
                 atoms=None,zmat=None,bonds=None,angles=None,dihedrals=None):
        self.name  = name
        self.backbone = backbone
        self.atoms = atoms
#        self.zmat = {}
        self.zmat = zmat
#        self.bonds = []
        self.bonds = bonds
#        self.angles = []
        self.angles = angles
#        self.dihedrals = []
        self.dihedrals = dihedrals
        self.rotate = rotate
        self.translate = translate
        
    def __str__(self):
        lines = []
        lines.append('mode template')
        lines.append('residue %s' % self.name)
        lines.append('info rotate %5.3f translate %5.3f' % (self.rotate,
                                                           self.translate))
        lines.append('backbone first %s middle %s last %s' % (
            self.backbone['first'],
            self.backbone['middle'],
            self.backbone['last']))
        atcount = 1
        found = True
        while(found):
            found = False
            for atName,atom in self.atoms.items():
                if atom['count'] == atcount:
                    found = True
                    string = ' '.join(atom['connec'])
                    par0 = atom['params']['middle']
                    par1 = par0
                    par2 = atom['params']['first']
                    par3 = par2
                    par4 = atom['params']['last']
                    par5 = par4
                    lines.append('atom %s %d %d %d %d %d %d %s ' % (atName,par0,par1,par2,par3,par4,par5,string) )
                    atcount = atcount + 1

#        for atName, atom in self.atoms.items():
#            string = ' '.join(atom['connec'])
#            par0 = atom['params']['middle']
#            par1 = par0
#            #What about par for first and last ??
#            lines.append('atom %s %d %d %s ' % (atName,par0,par1,string) )
        for atName, zmat in self.zmat.items():
            string = ' '.join(zmat)
            lines.append('zmat %s %s' % (atName,string) )
        for bond in self.bonds:
            string = ' '.join(bond)
            lines.append('bond %s ' % (string) )
        for angle in self.angles:
            string = ' '.join(angle)
            lines.append('angle %s ' % (string) )
        for dihedral in self.dihedrals:
            string = ' '.join(dihedral)
            lines.append('dihedral %s ' % (string) )            
        lines.append('\n')
        return '\n'.join(lines)
    
class Chain(object):
    def __init__(self,name='',bbAtoms=None,atoms=None,zmat=None,
                 bonds=None,angles=None,dihedrals=None):
        self.name = name
        self.bbAtoms = {}
        self.bbAtoms = bbAtoms
        self.atoms = {}
        self.atoms = atoms
        self.zmat = {}
        self.zmat = zmat
        self.bonds = []
        self.bonds = bonds
        self.angles = []
        self.angles = angles
        self.dihedrals = []
        self.dihedrals = dihedrals
        
    def __str__(self):
        lines = []
        lines.append('mode template')
        lines.append('chain %s' % self.name)
        
        lines.append('bbatom 1 %s %d %d ' % ('N',
                                             self.bbAtoms['N']['params'],
                                             self.bbAtoms['N']['params']))
        lines.append('bbatom 2 %s %d %d ' % ('CA',
                                             self.bbAtoms['CA']['params'],
                                             self.bbAtoms['CA']['params']))
        lines.append('bbatom 3 %s %d %d ' % ('C',
                                             self.bbAtoms['C']['params'],
                                             self.bbAtoms['C']['params']))
        lines.append('bbatom 4 %s %d %d ' % ('O',
                                             self.bbAtoms['O']['params'],
                                             self.bbAtoms['O']['params']))
        atcount = 1
        found = True
        while(found):
            found = False
            for atName,atom in self.atoms.items():
                if atom['count'] == atcount:
                    found = True
                    string = ' '.join(atom['connec'])
                    par0 = atom['params']
                    par1 = par0
                    lines.append('atom %s %d %d %s ' % (atName,par0,par1,
                                                        string) )
                    atcount = atcount + 1
##         for atName, atom in self.atoms.items():
##             string = ' '.join(atom['connec'])
##             par0 = atom['params']
##             par1 = par0
##             lines.append('atom %s %d %d %s ' % (atName,par0,par1,string) )
        for atName, zmat in self.zmat.items():
            string = ' '.join(zmat)
            lines.append('zmat %s %s' % (atName,string) )
        for bond in self.bonds:
            string = ' '.join(bond)
            lines.append('bond %s ' % (string) )
        for angle in self.angles:
            string = ' '.join(angle)
            lines.append('angle %s ' % (string) )
        for dihedral in self.dihedrals:
            string = ' '.join(dihedral)
            lines.append('dihedral %s ' % (string) )            
        lines.append('\n')            
        return '\n'.join(lines)
##         return """Chain %s
## bbAtoms %s
## atoms %s
## zmat %s
## bonds %s
## angles %s
## dihedrals %s
## """ % (self.name,self.bbAtoms,self.atoms,self.zmat,self.bonds,self.angles,
##        self.dihedrals)

def usage():
    print """USAGE : %s -c <file> -p <file> -h  

    -c / --connectivity= <file>
    template connectivity file
    
    -p / --parameters= <file>
    template parameters file

    -h / --help
    Prints this message
    
    """ % (sys.argv[0])
    sys.exit(-1)
    
def parse():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:p:h",
                                   ["connectivity=","parameters=","help"],)
    except getopt.GetoptError:
        usage()
    connecFile = None
    paramFile = None
    for o, a in opts:
        if o in ("-c", "--connectivity"):
            connecFile = a
        if o in ("-p", "--parameters"):
            paramFile = a
        if o in ("-h", "--help"):
            usage()
    if connecFile is None or paramFile is None:
        usage()
    return connecFile,paramFile

def readConnecFile(file=None):
    if file is None:
        raise "You must pass a filename!"
    stream = open(file,'r')
    buffer = stream.readlines()
    stream.close()
    # Process the contents of the file and generate template objects
    numTemplates = 0
    residues = {}
    chains = {}
    bbAtoms = {}
    atoms = {}
    zmat = {}
    backbone = {}
    bonds = []
    angles = []
    dihedrals = []
    for line in buffer:
        if line.startswith('mode'):
            if line.find('template') != -1:
                if numTemplates != 0:
                    if chain != '':
                        newChain = Chain(name=chain,
                                         bbAtoms=bbAtoms,
                                         atoms=atoms,
                                         zmat=zmat,
                                         bonds=bonds,
                                         angles=angles,
                                         dihedrals=dihedrals)
                        chains[chain] = newChain
                        #chains.append(newChain)
                    elif residue != '':
                        newResidue = Residue(name=residue,
                                             rotate=rotate,
                                             translate=translate,
                                             backbone=backbone,
                                             atoms=atoms,
                                             zmat=zmat,
                                             bonds=bonds,
                                             angles=angles,
                                             dihedrals=dihedrals)
                        residues[residue] = newResidue
                        #residues.append(newResidue)
                # Clear data
                numTemplates = numTemplates + 1
                chain = ''
                residue = ''
                bbAtoms = {}
                atoms = {}
                zmat = {}
                bonds = []
                angles = []
                atCount = 0
                dihedrals = []
                backbone = {}
                info = ''
        if line.startswith('chain'):
            elems = line.split()
            chain = elems[1].strip(' ')
        if line.startswith('residue'):
            elems = line.split()
            residue = elems[1]            
        if line.startswith('bbatom'):
            elems = line.split()
            bbAtoms[elems[2]] = {}
            bbAtoms[elems[2]]['params'] = {}
        if line.startswith('atom'):
            elems = line.split()
            atCount = atCount + 1
            atoms[elems[1]] = {}
            atoms[elems[1]]['connec'] = elems[2:]
            atoms[elems[1]]['params'] = {}
            atoms[elems[1]]['count'] = atCount
        if line.startswith('zmat'):
            elems = line.split()
            zmat[elems[1]]= elems[2:]            
        if line.startswith('bond'):
            elems = line.split()
            bonds.append(elems[1:])
        if line.startswith('angle'):
            elems = line.split()
            angles.append(elems[1:])
        if line.startswith('dihedral'):
            elems = line.split()
            dihedrals.append(elems[1:])
        if line.startswith('info'):
            elems = line.split()
            rotate = float(elems[2])
            translate = float(elems[4])
        if line.startswith('backbone'):
            elems = line.split()
            try:
                backbone[elems[1]] = elems[2]
                backbone[elems[3]] = elems[4]
                backbone[elems[5]] = elems[6]
                backbone[elems[7]] = elems[8]
            except IndexError:
                pass
    # Do the last one !
    newResidue = Residue(name=residue,
                         rotate=rotate,
                         translate=translate,
                         backbone=backbone,
                         atoms=atoms,
                         zmat=zmat,
                         bonds=bonds,
                         angles=angles,
                         dihedrals=dihedrals)
    residues[residue] = newResidue
    return chains, residues

def readParamFile(file=None):
    if file is None:
        raise "You must pass a file!"
    stream = open(file,'r')
    buffer = stream.readlines()
    stream.close()

    residues = {}
    inResidue = False
    inChain = False
    for line in buffer:
        if line.startswith('residue'):
            elems = line.split()
            name = elems[1]
            if name == 'HID':
                name = 'HIS'
            residues[name] = {}
            residues[name]['atoms'] = {}
            residues[name]['chain'] = {}
            
            inResidue = True
            inChain = False
            continue
        if line.startswith('chain'):
            elems = line.split()
            chain = elems[1]
            inChain = True
            residues[name]['chain'][chain] = {}
            continue
        if line.startswith('atom') and not inChain:
            elems = line.split()
            residues[name]['atoms'][elems[1]] = {}
            residues[name]['atoms'][elems[1]]['first'] = int(elems[3])
            residues[name]['atoms'][elems[1]]['middle'] = int(elems[5])
            residues[name]['atoms'][elems[1]]['last'] = int(elems[7])
        if line.startswith('atom') and inChain:
            elems = line.split()
            residues[name]['chain'][chain][elems[1]] = {}
            residues[name]['chain'][chain][elems[1]] = int(elems[2])
            residues[name]['chain'][chain][elems[1]] = int(elems[2])
            residues[name]['chain'][chain][elems[1]] = int(elems[2]) 
    return residues

def assignParameters(chains=None,residues=None,params=None):
    if chains is None or residues is None or params is None:
        raise "You must pass chains/residues/params"
    for name, residue in residues.items():
        #print name
        for atName, atom in residue.atoms.items():
            numbers = params[name]['atoms'][atName]
            residue.atoms[atName]['params'] = numbers
            #print atName
        for chainName, chain in residue.backbone.items():
            #print chain
            for atName, atom in chains[chain].bbAtoms.items():
                numbers = params[name]['chain'][chain][atName]
                chains[chain].bbAtoms[atName]['params'] = numbers
            for atName, atom in chains[chain].atoms.items():
                numbers = params[name]['chain'][chain][atName]
                chains[chain].atoms[atName]['params'] = numbers
        #sys.exit(-1)
#
# Actual script starts from here
#
if __name__ == '__main__':
    connecFile, paramFile = parse()
    chains, residues = readConnecFile(file=connecFile)
    # Build an amberchains dictionnary from chains
    # This is because the amber forcefield parameters can not be split
    # into backbone/sidechains the way OPLS works.
    amberChains = {}
    amberResidues = {}
    for name, residue in residues.items():
        #print name
        newBack = {}
        for position, chainName in residue.backbone.items():
            # Ignore single because no parameters are available in AMBER
            if position == 'single':
                continue
            #print position, chainName
            # Make new chain
            newChains = copy.deepcopy(chains[chainName])
            newChains.name = newChains.name+name
            amberChains[newChains.name] = newChains
            newBack[position] = newChains.name
        # Make new residue
        newResidue = copy.deepcopy(residue)
        newResidue.backbone = newBack
        amberResidues[newResidue.name] = newResidue
    # Now read parameter file for amber
    amberParams = readParamFile(file=paramFile)
    # Assign parameters to the residues/chains
    assignParameters(chains=amberChains,
                     residues=amberResidues,
                     params=amberParams)
    # Output everything
    for chain in amberChains.values():
        print chain
    for residue in amberResidues.values():
        print residue
        #sys.exit(-1)
