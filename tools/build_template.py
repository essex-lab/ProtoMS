# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Routines to build a ProtoMS template file

The module defines a two public function
build_template
make_zmat

Can be executed from the command line as a stand-alone program
"""

import sys
import os

import simulationobjects as sim

class PrepiAtom :
  """
  Class to encapsulate an atom in an Amber prepi file
  and its connectivity
  
  Attributes
  ----------
  atype : string
    the atom type
  bondidx : int
    the index of the atom this atom is bonded to in the prepi z-matrix
  bonds : list
    names of all atoms this atom is bonded to
  charge : float
  	the charge
  index : int
    the serial number in the prepi file
  is_on_loop : dictionary
    a flag for each atom in bonds indicating if the bond is part of a loop
  loop_closure : dictionary
    a flag for each atom in bonds indicating if the bond was part of a loop statement in the prepi file
  name : string
    the atom name
  traversed : dictionary
    a flag for each atom in bonds indicating if the bond has been traversed
  zmat : list
    atom names of the atoms defined to this atom in the z-matrix
  """
  def __init__(self) :
    self.charge = 0.0
    self.name = ""
    self.index = 0
    self.atype = ""
    self.bondidx = -1
    self.bonds = []
    self.loop_closure = {}
    self.traversed = {}
    self.zmat = []
    self.is_on_loop = {}
  def read(self,record,atomlist) :
    """
    Read a line from an Amber prepi-file
    """
    cols = record.strip().split()
    self.name = cols[1]
    self.index = int ( cols[0] )
    self.atype = cols[2]
    self.bondidx = int(cols[4])
    if self.bondidx >= 4 :
      self.add_bond(atomlist[self.bondidx-1])
    self.charge = float ( cols[10] )
  def add_bond(self,atomname,loop_closure=False) :
    """
    Add a bond to this atom
    """
    self.bonds.append(atomname)
    self.loop_closure[atomname] = loop_closure
    self.is_on_loop[atomname] = False
    self.traversed[atomname] = False
  def sort_bonds(self,metric) :
    """
    Sort the bonds based on some metric
    """
    self.bonds = sorted(self.bonds,key = lambda x : metric[x])[::-1]
  def next_bond(self,defined,update=True) :
    """
    Find an atom, bonded to this atom that 
    has not been defined and where the bond has not been traversed
    """
    next = None
    for bond in self.bonds :
      if not defined[bond] and not self.traversed[bond] and self.is_on_loop[bond]: #and not self.loop_closure[bond] :
        if update : self.traversed[bond] = True
        next = bond
        break
 
    if next == None :
      for bond in self.bonds :
        if not defined[bond] and not self.traversed[bond] and not self.is_on_loop[bond]: #and not self.loop_closure[bond] :
          if update : self.traversed[bond] = True
          next = bond
          break
    return next
  def backward_bond(self,exclude,defined) :
    """
    Find an atom, bonded to this atom that has
    already been defined but is not in the list of excluded atoms
    """
    back = None
    for bond in self.bonds :
      if bond not in exclude and defined[bond] :
        back = bond
        break
    return back
  def improper_dihedral(self,exclude,defined) :
    """
    Look for the definition of an improper dihedral:
    see if at least two bonds to this atom has been defined, 
    excluding the bond in the exclude list
    """
    a1 = None
    a2 = None
    for bond in self.bonds : 
      if bond in exclude : continue
      if defined[bond] :
        if a1 == None :
          a1 = bond
        elif a2 == None :
          a2 = bond
          return (a1,a2)
    return (a1,a2) 

def _find_cycles(atoms) :

  # Create a simple, unweighted, molecular graph
  graph = {}
  for atom in atoms :
    graph[atom] = {}
    for bond in atoms[atom].bonds : graph[atom][bond] = {}

  gnodes=set(graph)
  cycles=[]
  root = None
  while gnodes:  # loop over connected components
      if root is None:
          root=gnodes.pop()
      stack=[root]
      pred={root:root}
      used={root:set()}
      while stack:  # walk the spanning tree finding cycles
          z=stack.pop()  # use last-in so cycles easier to find
          zused=used[z]
          for nbr in graph[z]:
              if nbr not in used:   # new node
                  pred[nbr]=z
                  stack.append(nbr)
                  used[nbr]=set([z])
              elif nbr == z:        # self loops
                  cycles.append([z])
              elif nbr not in zused:# found a cycle
                  pn=used[nbr]
                  cycle=[nbr,z]
                  p=pred[z]
                  while p not in pn:
                      cycle.append(p)
                      p=pred[p]
                  cycle.append(p)
                  cycles.append(cycle)
                  used[nbr].add(z)
      gnodes-=set(pred)
      root=None
  return cycles

def _read_prepi(filename) :
  """
  Read an prepi-file from Antechamber
  """

  atomlist = "DUM DUM DUM".split()
  atoms = {}

  lines = open(filename,"r").readlines()
  i = 0
  while lines[i].find("CORRECT     OMIT") == -1 : i = i +1
  i = i + 5

  while len(lines[i]) > 4 :
    atom = PrepiAtom()
    atom.read(lines[i],atomlist)
    #graph.add_node(atom.name)
    atoms[atom.name] = atom
    atomlist.append(atom.name)
    if atom.bondidx >= 4 :
      atoms[atom.bonds[-1]].add_bond(atom.name)
    i = i + 1

  while lines[i].find("LOOP") == -1 : i = i + 1
  i = i +1
  while len(lines[i]) > 4 : 
    atom1,atom2 = lines[i].strip().split()
    atoms[atom1].add_bond(atom2,loop_closure=True)
    atoms[atom2].add_bond(atom1)
    i = i + 1
    
  # Find cycles in the atom graph
  cycles = _find_cycles(atoms)
  for cycle in cycles :
    for atom1,atom2 in zip(cycle[:-1],cycle[1:])  :
      atoms[atom1].is_on_loop[atom2] = True
      atoms[atom2].is_on_loop[atom1] = True
    atoms[cycle[0]].is_on_loop[cycle[-1]] = True
    atoms[cycle[-1]].is_on_loop[cycle[0]] = True
    
  return atomlist,atoms

#
# Calculate closeness centrality of all atoms in the molecular graph
#
def _closeness_centrality(graph) :
  # Return the length of the shortest path from node to all other nodes
  def shortest_path(graph,node) :
    seen={}
    level=0
    nextlevel={node:1}
    while nextlevel:
      thislevel=nextlevel
      nextlevel={} 
      for v in thislevel:
        if v not in seen:
          seen[v]=level
          nextlevel.update(graph[v])
      level=level+1
    return seen 
    
  c = {}
  for node in graph :
    sp = shortest_path(graph,node) 
    spsum = sum(sp.values())
    # I don't believe the normalization is necessary since every atom can be reached from evert other atom
    #norm = (len(sp)-1.0) / (len(graph) - 1.0)
    c[node] = (len(sp) - 1.0) / spsum
  return c

#
# Compute the closeness of all atoms, sort
# their bonds on closeness and return of sorted list of atoms
#
def _compute_closeness(atoms,verbose=False) :

  # Create a simple, unweighted, molecular graph
  graph = {}
  for atom in atoms :
    graph[atom] = {}
    for bond in atoms[atom].bonds : graph[atom][bond] = {}

  # Calculate the closeness of all atoms
  closeness = _closeness_centrality(graph)

  # Create a sorted list based on the closeness of all atoms
  atomlist = sorted(closeness,key=closeness.get)[::-1]

  # Sort the bonds for each atom by the closeness and set all atoms to undefined
  for atom in atomlist :
    atoms[atom].sort_bonds(closeness)
    if verbose : print "%s %.3f"%(atom,closeness[atom])
  if verbose : print ""
    
  return atomlist

def make_zmat(prepifile):
  """ Make a ProtoMS z-matrix for a solute
   
  make_zmat(prepifile)
  
  Parameters
  ----------   
  prepifile : string
    the filename of the Amber prepi file (prep-file with z-matrix)

  Returns
  -------
  atoms : dictionary
    created PrepiAtom objects
  all_zmat : list
    a list of atoms with their z-matrix atoms
  """
    
  # ----------------
  # Helper routines
  # ----------------

  def define_atom(current,previous,atoms,defined,all_zmat,verbose=False) :
    """
    Define a new atom in the z-matrix
    """
  # Make a proper dihedral for the current atom
    def make_proper() :
      # If we a dealing with the first three atoms, just include less and less dummies
      if len(all_zmat) <= 2 : 
        atoms[current].zmat = atoms[previous].zmat[:2]
        atoms[current].zmat.insert(0,previous)
      # if not, we traverse backwards on the molecular graph, looking for two atoms to
      # form an angle and dihedral
      else :
        a2 = atoms[previous].backward_bond([current,previous],defined)
        a3 = atoms[a2].backward_bond([previous,a2],defined)
        atoms[current].zmat = ("%s %s %s"%(previous,a2,a3)).split()
      
    # Make an improper dihedral for the current atom
    def make_improper(a2,a3) :
      # This two lines were the previous default definition, now they are supplied as arguments
      #a2 = atoms[previous].bonds[0]
      #a3 = atoms[previous].bonds[1]
      atoms[current].zmat = ("%s %s %s"%(previous,a2,a3)).split()
          
    defined[current] = True
    if previous != None :
    
      if len(all_zmat) <= 2 : 
        make_proper()
      else :
  
        # This was the previous rule:
        # If the current atom is the second or more bonded atom to the previous,
        # let's define it with an improper so that it moves with the first atom
        # bonded to the previous atom
        #if atoms[previous].bonds.index(current) > 1 :
    
        # Now, we let the previous atom make the decision:
        # simply return two previously defined atoms, or None
        a2,a3 = atoms[previous].improper_dihedral(current,defined)
        if a2 != None and a3 != None :
          make_improper(a2,a3)
        else :
          make_proper()
 
    # This will only happens for the very first atom
    else :
      atoms[current].zmat = "DM3 DM2 DM1".split()
    if verbose : print " ",current,
    all_zmat.append(current+" "+" ".join(atoms[current].zmat))

  def traverse_graph(atomlist,atoms,verbose=False) :
    """
    Traverse the molecule graph based on the closeness of the atoms
    and define the z-matrix
    """
    defined = {}
    for atom in atoms : defined[atom] = False
    all_zmat = []
  
    if verbose : print "Traversal of the molecular graph:"
  
    # Start with the most central atom
    branch_atom = atomlist[0]
    define_atom(branch_atom,None,atoms,defined,all_zmat,verbose)
  
    while True :
      # Traverse a branch from an atom with at least two bonds
      previous = branch_atom
      next = atoms[branch_atom].next_bond(defined)
      atoms[next].traversed[branch_atom] = True
      while next != None : 
        define_atom(next,previous,atoms,defined,all_zmat,verbose)
        previous = next
        next = atoms[previous].next_bond(defined)
        if next != None : 
          atoms[next].traversed[previous] = True
      if verbose : print ".",
    
      # Check if we have more branches to traverse
      if atoms[branch_atom].next_bond(defined,update=False) == None :
        if verbose : print " : ",
        # Tries to find a new atom to branch off from
        found = False
        for atom  in atomlist :
          if not defined[atom] :
            branch_atom = atom
            define_atom(branch_atom,atoms[atom].bonds[0],atoms,defined,all_zmat,verbose)
            # If the found atom has more than one bonds, we can traverse its branches
            if len(atoms[branch_atom].bonds) > 1 :
              found = True
              break
            else :
              if verbose : print "; ",
        # If we could not find any more atoms to branch off from we are done!      
        if not found : break 
    if verbose : print ""
  
    return all_zmat
    
  # Parse the prepi-file into a list of atom names and a dictionary of Atom objects
  atomnames,atoms = _read_prepi(prepifile)	# MOD
  h,t = os.path.splitext(prepifile)	

  # Compute closeness of all atoms and sort their bonds based on this  
  atomlist_closeness = _compute_closeness(atoms,verbose=False) 		# MOD
  # Traverse the molecular graph and define the z-matrix
  all_zmat = traverse_graph(atomlist_closeness,atoms,verbose=False)	# MOD
    
  #with open ( '%s.zmat' % h, 'w' ) as f:
  #  f.write ( '\n'.join ( all_zmat ) )
        
  return atoms,all_zmat

def _readfrcmod(filename):
    """ Read an Amber frcmod file from disc
   
    readfrcmod(filename)
  
    Parameters
    ----------   
    filename - the filename of the Amber frcmod file
    
    Returns
    -------
    a list of bond, angle and dihedral parameters
    """
    with open ( filename ) as f:
        #Find start of bond region
        while not f.next().startswith ( 'BOND' ):
            pass

        bonds = []
        for line in f:
            cols = line.split ()
            if not cols:
                break
            bonds.append ( [ cols[0].split ( '-' ) ] + map ( float, cols[1:3] ) )
        f.next()

        angles = []
        for line in f:
            cols = line.split ()
            if not cols:
                break
            angles.append ( [ cols[0].split ( '-' ) ] + map ( float, cols[1:3] ) )

        f.next()

        dihedrals = []
        for line in f:
            cols = line.split ()
            if not cols:
                break
            dihedrals.append ( [ cols[0].split ( '-' ) ] + map ( float, cols[1:5] ) )

        f.next()
            
    return bonds, angles, dihedrals



def build_template ( temfile, prepifile, translate=0.25, rotate=5, zmatfile=None, frcmodfile=None, resname="UNK" ) :
    """ Build a ProtoMS template file
   
    build_template ( temfile, prepifile, zmatfile=None, frcmodfile=None, resname="UNK" ) 
  
    Parameters
    ----------   
    temfile - the filename to save the template file to
    prepifile - the filename of the Amber prepi file (prep-file with z-matrix)
    zmatfile - the filename of a zmat, if None it is created
    frcmodfile - the filename of an Amber frcmod file with additional parameters
    resname - the name of solute

    Returns
    -------
    None
    """
    if zmatfile is None :
        atoms,zmat = make_zmat(prepifile)
    else :
        with open ( zmatfile ) as f:
            zmat = f.readlines()
        atoms = _read_prepi ( prepifile )[1]

    centralities = _compute_closeness ( atoms, verbose = False )

    if frcmodfile is None :
      frcbonds = frcangles = frcdihedrals = None
    else :  
      frcbonds, frcangles, frcdihedrals = _readfrcmod ( frcmodfile )

    PROTOMSHOME = os.getenv('PROTOMSHOME')

    angle_params = sim.parameter_set ( "%s/parameter/gaff14.ff" % PROTOMSHOME, 'angle' )
    dihedral_params = sim.parameter_set ( "%s/parameter/gaff14.ff" % PROTOMSHOME, 'dihedral' )

    with open ( "%s/parameter/gaff.types" % PROTOMSHOME ) as f:
        at_params = [ line.split() for line in f ]

    kBT = 0.0019872041 * 300 #Boltzmann constant at 300 kelvin in kcal/mol
    move_scale = 0.5

    dummynames = ["DM3", "DM2", "DM1"]
    aromatic_types = "ca cp cq ce cf cc cd nb ne nf pb pe pf px py sx sy".split()
    template = sim.TemplateFile()
    template.templates.append(sim.MolTemplate())
    template.templates[0].atomclass = sim.TemplateSoluteAtom
    moltem = template.templates[0]
    moltem.type = "solute"
    moltem.name = resname
    moltem.translate = translate
    moltem.rotate = rotate
    
    def get_resname ( s ):
        if s in dummynames:
            return 'DUM'
        else:
            return resname

    out = ''
    if frcbonds is not None :
        out += """mode bond
# U(r) = k(r-r0)**2
#parameter k(kcal mol-1 A-2) r0(A) comment\n"""
        for i, bond in enumerate ( frcbonds, 4500 ):
            #out += "par %4d %.3f %.3f\n" % ( i, bond[1], bond[2] )
            template.bondparams.append(sim.ForceFieldParameter(record="par %4d %.3f %.3f\n" % ( i, bond[1], bond[2] )))
        if not frcbonds:
            out += '\n'
        out += "#atm atm1 atm2 parameter\n"
        for i, bond in enumerate ( frcbonds, 4500 ):
            out += "atm %s %s %4d\n" % ( bond[0][0], bond[0][1], i )
            template.bondatoms.append(sim.AtomSet(record="atm %s %s %4d\n" % ( bond[0][0], bond[0][1], i )))
        if not frcbonds:
            out += '\n'

    if frcangles is not None :
        out += """mode angle
# U(theta) = k(theta-theta0)**2
#parameter k(kcal mol-1 deg-2) theta0(deg) comment\n"""
        for i, angle in enumerate ( frcangles, 4700 ):
            out += "par %4d %.3f %.3f\n" % ( i, angle[1], angle[2] )
            template.angleparams.append(sim.ForceFieldParameter(record="par %4d %.3f %.3f\n" % ( i, angle[1], angle[2] )))
        if not frcangles:
            out += '\n'
        out += "#atm atm1 atm2 parameter\n"
        for i, angle in enumerate ( frcangles, 4700 ):
            out += "atm %s %s %s %4d\n" % ( angle[0][0], angle[0][1], angle[0][2], i )
            template.angleatoms.append(sim.AtomSet(record="atm %s %s %s %4d\n" % ( angle[0][0], angle[0][1], angle[0][2], i )))
        if not frcangles:
            out += '\n'

    if frcdihedrals is not None :
        out += """mode dihedral
# U(phi) = k1*( 1.0 + k2*(cos[k3*phi + k4]) )
#term k1(kcal mol-1) k2 k3 k4(deg) #comment\n"""
        fmt = ( 4600, 0.0, 0.0, 0.0, 0.0 )
        out += "term  %4d   %.3f   %.3f   %.3f   %.3f\n" % fmt
        for i, di in enumerate ( frcdihedrals, 4601 ):
            fmt = ( i, di[2], di[1], di[4], di[3] )
            out += "term  %4d   %.3f   %.3f   %.3f   %.3f\n" % fmt
            template.dihedralterms.append(sim.ForceFieldParameter(record="term  %4d   %.3f   %.3f   %.3f   %.3f\n" % fmt))
        out += "#par  term1  term2 etc..  #comment\n"

        #Loop through unique dihedral terms
        done = []
        for i, di in enumerate ( frcdihedrals, 4600 ):
            if di[0] not in done:
                #Get all terms with matching atom groups from dihedral list
                terms = [ ( j, term ) 
                          for j, term in enumerate ( frcdihedrals, 4601 ) 
                          if term[0] == di[0] ]
                #Start par entry
                out += "par  %4d" % i
                record = "par  %4d" % i
                #Loop through all 5 terms of full dihedral
                for j in xrange ( 5 ):
                    #Pull out term with values of k3 == j
                    ref = [ term for term in terms if int ( term[1][4] ) == j ]
                    if ref:
                        out += "   %4d" % ref[0][0]
                        record += "   %4d" % ref[0][0]
                    else:
                        #If no k3 == j give default zero parameters
                        out += "   %4d" % 4600
                        record += "   %4d" % 4600
                out += '\n'
                template.dihedralparams.append(sim.ForceFieldParameter(record=record))
            done += [ di[0] ]
        if not frcdihedrals:
            out += '\n'
        out += "#atm atm1 atm2 atm3 atm4 parameter #comment\n"

        done = []
        for i, di in enumerate ( frcdihedrals, 4600 ):
            if di[0] not in done:
                fmt = ( di[0][0], di[0][1], di[0][2], di[0][3], i )
                out += "atm %4s %4s %4s %4s %4d\n" % fmt
                template.dihedralatoms.append(sim.AtomSet(record="atm %4s %4s %4s %4s %4d\n" % fmt))
            done += [ di[0] ]
        if not frcdihedrals:
            out += '\n'

    out += """mode clj
#parameter atm proton-num charge(|e|) sigma(A) epsilon(kcal mol-1)\n"""
    # Here we need to iterate over the atom order in the z-matrix, so this is modified
    for i, atom in enumerate ( [atoms[zatms.split()[0]] for zatms in zmat], 3000 ):
        params = [ j for j in at_params if j[0] == atom.atype ][0]
        fmt = ( i, atom.atype, int ( params[6] ), 
                atom.charge, float ( params[2] ), float ( params[4] ) )
        out += "par  %4d %4s   %02d %10.5f %10.5f %10.5f\n" % fmt
        template.cljparams.append(sim.ForceFieldParameter(record="par  %4d %4s   %02d %10.5f %10.5f %10.5f\n" % fmt))
    # out += '\n'

    #out += zmat
    out += """mode template
solute %s
info translate %f rotate %f\n""" % ( resname, translate, rotate )

    # Print out the atoms
    for i, line in enumerate ( zmat, 3000 ):
        atms = line.split()
        fmt = ( atms[0], resname, i, i, 
                atms[1], get_resname ( atms[1] ), 
                atms[2], get_resname ( atms[2] ), 
                atms[3], get_resname ( atms[3] ) )
        out += "atom  %4s %4s %4i %4i %4s %4s %4s %4s %4s %4s\n" % fmt
        moltem.atoms.append(sim.TemplateSoluteAtom(record="atom  %4s %4s %4i %4i %4s %4s %4s %4s %4s %4s\n" % fmt))
        
    # Print out the bonds
    taken_bonds = []
    for line in zmat[1:]:
        atom = line.split()[:2]
        taken_bonds.append((atom[0],atom[1]))
        taken_bonds.append((atom[1],atom[0]))
        out += "bond %4s %3s %4s %3s\n" % (atom[0],resname,atom[1],resname)
        moltem.connectivity.append(sim.TemplateConnectivity(record="bond %4s %3s %4s %3s\n" % (atom[0],resname,atom[1],resname)))
        
    # Also need to do bonds that close loops
    for atom1 in atoms :
      for atom2 in atoms[atom1].bonds :
        if (atom1,atom2) not in taken_bonds :
          taken_bonds.append((atom1,atom2))
          taken_bonds.append((atom2,atom1))
          out += "bond %4s %3s %4s %3s\n" % (atom1,resname,atom2,resname)
          moltem.connectivity.append(sim.TemplateConnectivity(record="bond %4s %3s %4s %3s\n" % (atom1,resname,atom2,resname)))
          
# THIS CODE IF BROKEN, replaced with the one above
#    for atom in atoms.itervalues():
#        for i in atom.loop_closure:
#            if atom.loop_closure[i] == True:
#                out += "bond %4s %3s %4s %3s\n" % (atom.name,resname,i,resname)

    # Print out the angles:
    for line in zmat[2:]:
        angle = line.split()[:3]
        atypes = [ atoms[i].atype for i in angle ]
        try:
            k = angle_params.get_params ( atypes ).k
        except IndexError:
            #parameter not in gaff.ff try frcmod params
            try:
                k = [ i for i in frcangles
                      if atypes in ( i[0], i[0][::-1] ) ][0][1]
            except IndexError:
                print "Unable to find angle parameters for %s-%s-%s" % tuple ( atypes )
                print "For now the flexibility of this angle will be set to zero."
                print "To correct this consider manually adding a term to the frcmod file\n"
                k = 10**10 #make k huge so max_move is zero
            except TypeError :
                k = 10*10

        max_move = move_scale * 2 * ( 1 / ( kBT * k ) )**0.5
        if False in [ i in aromatic_types for i in atypes ]:
            fmt = ( angle[0], resname, 
                    angle[1], resname, 
                    angle[2], resname, max_move ) 
            out += "angle %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt
            moltem.connectivity.append(sim.TemplateConnectivity(record="angle %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt))
            
    min_flex = 2.0
    max_flex = 10.0
    for line in zmat[3:]: #First 3 lines contain dummies
        # The first 4 atoms in the z matrix corresponds to the dihedrals.
        dihedral = line.split()[:4]
        atypes = [ atoms[i].atype for i in dihedral ]
        # Check matching parameters exist in gaff.ff
        missing = False
        try:
            dihedral_params.get_params ( atypes )
        except IndexError:
            # If not in gaff.ff then check frcmod
            try:
                [ i for i in frcdihedrals
                  if atypes in ( i[0], i[0][::-1] ) ]
            except IndexError:
                print "Unable to find dihedral parameters for %s-%s-%s-%s" % tuple ( atypes )
                print "For now the flexibility of this angle will be set to zero."
                print "To correct this consider manually adding a term to the frcmod file\n"
                k = 10**10 #make k huge so max_move is zero
                missing = True
            
        # Ensure improper dihedrals are removed by checking the last 
        # atom in the dihedral is bonded to the previous atom in the list.
        if not dihedral[2] in atoms[dihedral[3]].bonds:
            continue

        # If central bond is part of a loop
        if atoms[dihedral[1]].is_on_loop[dihedral[2]]:
            # Ensure aromatic dihedrals are excluded
            if atypes[1] in aromatic_types and atypes[2] in aromatic_types:
                continue
            fmt = ( dihedral[0], resname, 
                    dihedral[1], resname, 
                    dihedral[2], resname, 
                    dihedral[3], resname, min_flex )
            out += "dihedral %4s %4s %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt
            moltem.connectivity.append(sim.TemplateConnectivity(record="dihedral %4s %4s %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt))
        else:
            # Not part of loop so base flexibility on centrality
            av_rank = ( centralities.index ( dihedral[1] ) + centralities.index ( dihedral[2] ) ) / float ( 2 * len ( centralities ) )
            flex = ( av_rank ) * ( max_flex - min_flex ) + min_flex
            # print centralities.index ( dihedral[1] ) + centralities.index ( dihedral[2] ), av_rank
            if missing:
                flex = 0.0
            fmt = ( dihedral[0], resname, 
                    dihedral[1], resname, 
                    dihedral[2], resname, 
                    dihedral[3], resname,flex)
            out += "dihedral %4s %4s %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt
            moltem.connectivity.append(sim.TemplateConnectivity(record="dihedral %4s %4s %4s %4s %4s %4s %4s %4s flex %.3f\n" % fmt))
            
    #with open ( temfile, 'w' ) as f:
    #    f.write ( out )
    #template.write(temfile)
    return template

if __name__ == '__main__':

    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to build a ProtoMS template file")
    parser.add_argument('-p','--prepi',help="the name of the leap prepi-file")
    parser.add_argument('-o','--out',help="the name of the template file",default="lig.tem")
    parser.add_argument('-z','--zmat',help="the name of the zmatrix-file, if it exists")
    parser.add_argument('-f','--frcmod',help="the name of the frcmod-file, if it exists")
    parser.add_argument('-n','--name',help="the name of the solute",default='UNK')
    parser.add_argument('-t','--translate',help="maxmium size for translation moves in Angstroms", default=0.25, type=float)
    parser.add_argument('-r','--rotate',help="maxmium size for rotation moves in degrees", default=5.0, type=float)
    args = parser.parse_args()
    
    tem = build_template ( temfile=args.out, prepifile=args.prepi, zmatfile=args.zmat, frcmodfile=args.frcmod, resname=args.name, translate=args.translate, rotate=args.rotate ) 
    