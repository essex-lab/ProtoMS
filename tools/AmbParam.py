"""Library of classes for dealing with Amber parameter files."""

import os, sys

class AmbParameterAtm:
    def __init__ ( self, args ):
        self.at = args[0]
        self.mass = int ( args[1] )

    def add ( self, args ):
        """Routine to add additional parameters for atoms at the end of the param file"""
        self.sigma = args[1]
        self.epsilon = args[2]
        
class AmbParameterBnd:
    def __init__ ( self, args ):
        self.at1 = args[0]
        self.at2 = args[1]
        self.k = args[2]
        self.b0 = args[3]
          
class AmbParameterAng:
    def __init__ ( self, args ):
        self.at1 = args[0]
        self.at2 = args[1]
        self.at3 = args[2]
        self.k = args[3]
        self.b0 = args[4]

class AmbTermDih:
    def __init__ ( self, args ):
        self.k1, self.k2, self.k3, self.k4 = args

        
class AmbParameterDih:
    def __init__ ( self, args ):
        self.at1 = args[0]
        self.at2 = args[1]
        self.at3 = args[2]
        self.at4 = args[3]
        self.terms = {}
        self.terms[abs(args[7])] = AmbTermDih ( args[4:] )

    def add ( self, args ):
        """Add term to this dihedral"""
        self.terms[abs(args[7])] = AmbTermDih ( args[4:] )


class AmbRes:
    """Class to store Amber residue information from .in parameter files"""
    def __init__ ( self ):
        self.ats = {}
        self.charges = {}

    def add_atom ( self, name, at, charge ):
        self.ats[name] = at
        self.charges[name] = charge

class clj:
    """Combination of Coloumb, Lennard-Jones and atom type information
    as required for a clj entry in a ProtoMS parameter file."""
    def __init__ ( self, id, name, type, chg, at ):
        # print name, chg, at.at
        self.id = id
        self.name = name
        self.type = type
        self.chg = chg
        self.sigma = at.sigma * 1.7818
        self.epsilon = at.epsilon
        
class cljs:
    """A class to handle collections of clj obeects with look up routines etc."""
  
    def __init__ ( self ):

        self.reses = {}
        self.curr_id = 1
        self.dum_at = AmbParameterAtm ( ( 'DU', 0 ) )
        self.dum_at.add ( ( 0, 0, 0 ) )
        self.dum = clj ( 100, 'DU', 'DU', 0, self.dum_at )
        self.pars = [ self.dum ]
        
    def dup ( self, par ):
        for i in self.pars:
            if ( par.type == i.type and par.chg == i.chg and
                 par.sigma == i.sigma and par.epsilon == i.epsilon ):
                return True, i
              
        return False, par

    def __iter__ ( self ):
        return iter ( sorted ( self.pars, key = lambda x: x.id ) )
        
    def add ( self, at, chg, name, type, res ):
        par = clj ( self.curr_id, name, type, chg, at )

        novel, par = self.dup ( par )
        if not novel:
            self.pars.append ( par )
            self.curr_id += 1
            #special fix as id 100 is a set of dummy parameters
            #by ProtoMS convention
            if self.curr_id == 100:
                self.curr_id += 1
            
        if res not in self.reses:
            self.reses[res] = {}

        self.reses[res][name] = par

    def cljbyres ( self, res, name ):
        """ """
        return self.reses[res][name]

class AmbParameterSet:
    """Class to handle AMBER parameter sets for automated production of ProtoMS parameter files.
    Routines for this class should be used in the following way.

    Amber parameter files should be parsed using the read routines.
    The fix_diffs routine should then be run to account for differences in AMBER and ProtoMS templates.
    The write_protoms_ff and write_protoms_residues will then produce a pair of valid ProtoMS parameter files.
 
    Amber parameter files must be read in the following order:
    .dat
    .frcmod
    .in
    """

  
    def __init__ ( self ):

        self.atms = {}
        self.angs = {}
        self.bnds = {}
        self.dihs = {}
        self.cljs = cljs ()

        self.atm_fmt = str, float, float
        self.bnd_fmt = str, str, float, float
        self.ang_fmt = str, str, str, float, float
        self.dih_fmt = str, str, str, str, float, float, float, float
        self.dih_fmt2 = str, str, str, str, float, float, float        

        self.Ns = 'N', 'NA', 'N2', 'N*', 'NC', 'NB', 'NT', 'NY'
        self.Cs = 'C*', 'CA', 'CB', 'CC', 'CD', 'CK', 'CM', 'CN', 'CQ', 'CR', 'CV', 'CW', 'CY', 'C'

    def read_in ( self, filename, term = None ):
        """Read in data from an Amber .in parameter file, as used with leap.
           Must be read in after .dat files"""

        if term == None or term == 'center':
            suffix = ''
        elif term == 'cterm':
            suffix = 'ct'
        elif term == 'nterm':
            suffix = 'nt'
        else:
            raise Exception ( 'Invalid value given as term argument. Must be None, center, cterm or nterm' )
        with open ( filename ) as f:
            for line in f:
                cols = line.strip().split()
                if len ( cols ) == 3 and cols[1] == 'INT':
                    res = cols[0] + suffix
                    # self.reses[res] = AmbRes ()
                if len ( cols ) == 11:
                    name = cols[1]
                    at_type = cols[2]
                    if at_type in ( 'DU', 'IP', 'IM' ):
                        continue
                    if at_type in self.Ns:
                        at2 = 'N'
                    elif at_type in self.Cs:
                        at2 = 'C'
                    else:
                        at2 = at_type
                    # clj = ( cols[1], self.atms[at2], float ( cols[10] ), at )
                    self.cljs.add ( self.atms[at2], float ( cols[10] ), name, at_type, res )
                  
    def read_dat ( self, filename ):
        """Routine to read in data from an Amber .dat parameter file, as used with leap."""

        added = []
        with open ( filename ) as f:
            for line in f:
                if line == '\n':
                    continue
                line = line.replace ( '-', ' ' )
                cols = line.strip().split()

                #remove dashes so we can interpret lines easily
                line = line.replace ( '-', ' ' )
                cols = line.strip().split()

                #logic for parsing parameter files adapted from leap src code
                #
                #if the parameter line can be mapped to the appropriate format combination
                #without throwing an exception and returns enough arguments then
                #its the right parameter type
                #Atom Type
                try:
                    args = tuple ( i ( j ) for i, j in zip ( self.atm_fmt, cols ) )
                except ValueError:
                    pass
                else:
                    if len ( args ) < 3:
                        continue
                    #some extra logic required to deal with recuring atom entries
                    if args[0] not in added:
                        self.atms[args[0]] = AmbParameterAtm ( args )
                    else:
                        self.atms[args[0]].add ( args )
                    added.append ( args[0] )
                    continue

                #Bond
                try:
                    args = tuple ( i ( j ) for i, j in zip ( self.bnd_fmt, cols ) )
                except ValueError:
                    pass
                else:
                    if len ( args ) < 4:
                        continue
                    self.bnds[args[:2]] = AmbParameterBnd ( args )
                    continue

                #Angle
                try:
                    args = tuple ( i ( j ) for i, j in zip ( self.ang_fmt, cols ) )
                except ValueError:
                    pass
                else:
                    if len ( args ) < 5:
                        continue
                    self.angs[args[:3]] = AmbParameterAng ( args )
                    continue

                #Dihedral
                try:
                    args = tuple ( i ( j ) for i, j in zip ( self.dih_fmt, cols ) )
                except ValueError:
                    pass
                else:
                    if len ( args ) < 8:
                        continue
                    #some extra logic required to deal with multiple dihedral terms
                    if args[:4] not in added:
                        self.dihs[args[:4]] = AmbParameterDih ( args )
                    else:
                        self.dihs[args[:4]].add ( args )
                    # print args
                    added.append ( args[:4] )
                    continue

                
                  
                # #Some dihedrals miss a column for some reason - they are impropers
                # try:
                #     args = tuple ( i ( j ) for i, j in zip ( self.dih_fmt2, cols ) )
                # except ValueError:
                #     pass
                # else:
                #     if len ( args ) < 7:
                #         continue
                #     #some extra logic required to deal with multiple dihedral terms
                #     args = args[:4] + (1,) + args[4:]
                #     if args[:4] not in self.dihs:
                #         self.dihs[args[:4]] = AmbParameterDih ( args )
                #     else:
                #         self.dihs[args[:4]].add ( args )
                #     # print args
                #     continue

    def read_frcmod ( self, filename ):
      """Routine to read in data from an Amber .frcmod parameter file, as used with leap.
      Frcmods should be read in after the relevant .dat files to prevent overwriting of parameter sets.

      Frcmod files follow the same formatting rules as .dat so this function is simply a wrapper to read_dat
      provided for semantic reasons."""

      self.read_dat ( filename )

      
    def fix_diffs ( self ):
        """Apply special fixes to account for fundamental differences in amber and ProtoMS parameter files """

        #remove entries for TIP3P, these are handled separately in solvents.ff
        for i in self.bnds.keys():
            if 'OW' in i:
                self.bnds.pop ( i )

        #create new entries for n and c terminal GLH and ASH
        glh = self.cljs.reses['GLH']
        self.cljs.reses['GLHnt'] = { i : glh[i] for i in glh }
        self.cljs.reses['GLHct'] = { i : glh[i] for i in glh }

        ash = self.cljs.reses['ASH']
        self.cljs.reses['ASHnt'] = { i : ash[i] for i in ash }
        self.cljs.reses['ASHct'] = { i : ash[i] for i in ash }

        #Add/remove additional required atoms
        self.cljs.reses['GLHct']['OXT'] = self.cljs.reses['ALA']['O']
        for i in xrange ( 1, 4 ):
            self.cljs.reses['GLHnt']['H%d' % i] = glh['H']
        self.cljs.reses['GLHnt'].pop ( 'H' )

        self.cljs.reses['ASHct']['OXT'] = self.cljs.reses['ALA']['O']
        for i in xrange ( 1, 4 ):
            self.cljs.reses['ASHnt']['H%d' % i] = ash['H']
        self.cljs.reses['ASHnt'].pop ( 'H' )
        
        #special fix for HID due to diff in amber and ProtoMS naming
        for suff in '', 'nt', 'ct':
            self.cljs.reses['HIS' + suff] = self.cljs.reses.pop ( 'HID' + suff )

        #special fix for PRO due to missing hydrogen in backbone for nterm
        pro = self.cljs.reses['PROnt']
        pro['H1']= pro.pop ( 'H2' )
        pro['H2']= pro.pop ( 'H3' )

        #special fix for GLU to add dummy
        for suff in '', 'nt', 'ct':
            self.cljs.reses['GLU'+suff]['HE2'] = self.cljs.dum

        #special fix for HIS to add dummy
        for suff in '', 'nt', 'ct':
            self.cljs.reses['HIS'+suff]['HNE'] = self.cljs.dum

        #special fix for lysine
        lyn = self.cljs.reses['LYN']
        lyn['HZ1'] = self.cljs.dum
        self.cljs.reses['LYNnt'] = { i : self.cljs.dum for i in lyn }
        self.cljs.reses['LYNct'] = { i : self.cljs.dum for i in lyn }

        #special fix for cystiene, histidine and asparate
        for suff in '', 'nt', 'ct':
            self.cljs.reses['CYX'+suff]['HG'] = self.cljs.dum
            self.cljs.reses['HIE'+suff]['HD1'] = self.cljs.dum
            self.cljs.reses['ASP'+suff]['HD2'] = self.cljs.dum
        
            
    def write_protoms_ff ( self, fname ):

        s = """mode info
ljcombine arithmetic
scl14coul  0.833333
scl14lj    0.500\n"""
        
        s += '\nmode bond\n'
        for i, j in enumerate ( sorted ( self.bnds ), 1 ):
            bnd = self.bnds[j]
            s += 'par  %4d%8.1f%15.10f  #%s-%s\n' % ( i, bnd.k, bnd.b0, bnd.at1, bnd.at2 )

        for i, j in enumerate ( sorted ( self.bnds ), 1 ):
            bnd = self.bnds[j]
            s += 'atm  %4s  %4s  %4d #\n' % ( bnd.at1, bnd.at2, i )

        s += 'mode angle\n'
        for i, j in enumerate ( sorted ( self.angs ), 1 ):
            ang = self.angs[j]
            s += 'par  %4d%8.1f%8.3f  #%s-%s-%s\n' % ( i, ang.k, ang.b0, ang.at1, ang.at2, ang.at3 )

        for i, j in enumerate ( sorted ( self.angs ), 1 ):
            ang = self.angs[j]
            s += 'atm  %4s  %4s  %4s  %4d\n' % ( ang.at1, ang.at2, ang.at3, i )


        s += 'mode dihedral\n'
        c = 1
        for i in sorted ( self.dihs ):
            dih = self.dihs[i]
            s += 'term  %4d  %15.8f%8.3f%8.3f%8.3f\n' % ( c, 0, 0, 0, 0 )
            c += 1
            for j in xrange ( 1, 5 ):
                try:
                    term = dih.terms[ float ( j ) ]
                    s += 'term  %4d  %15.8f%8.3f%8.3f%8.3f\n' % ( c, term.k2 / term.k1, 1, j, term.k3 )
                except KeyError:
                    s += 'term  %4d  %15.8f%8.3f%8.3f%8.3f\n' % ( c, 0, 1, j, 0 )
                c += 1
                
        for i, j in enumerate ( sorted ( self.dihs ), 1 ):
            dih = self.dihs[j]
            terms = tuple ( range ( ( i - 1 ) * 5 + 1, ( i * 5 + 1  ) ) )
            s += 'par  %4d  %4d  %4d  %4d  %4d  %4d  # %s-%s-%s-%s\n' % ( ( i, ) + terms + ( dih.at1, dih.at2, dih.at3, dih.at4 ) )
            
        for i, j in enumerate ( sorted ( self.dihs ), 1 ):
            dih = self.dihs[j]
            s += 'atm  %4s  %4s  %4s  %4s  %4d  #\n' % ( dih.at1, dih.at2, dih.at3, dih.at4, i )

        prot_no = { 'H' : 1, 'C' : 6, 'N' : 7, 'O' : 8, 'S' : 16, 'D' : 0 }
        s += 'mode clj\n'
        for i, j in enumerate ( self.cljs, 1 ):
            element = j.type[0]
            try:
                float ( element )
                element = j.type[1]
            except ValueError:
              pass

            s += 'par  %4d  %4s   %02d  %8.4f  %8.4f  %8.4f  #\n' % ( ( i, j.type, prot_no[element], j.chg, j.sigma, j.epsilon ) )

        with open ( fname, 'w' ) as f:
            f.write ( s )

    def write_protoms_residues ( self, fname ):
        home = os.environ['PROTOMSHOME']
        sys.path.append ( os.path.join ( home, 'tools' ) )
        from convertatomnames import read_convfile
        conv = read_convfile ( os.path.join ( home, 'data', 'atomnamesmap2.dat' ),
                               inmode = 'protoms', outmode = 'amber' )
        conv['GLH'] = conv['GLU']
        
        s = ''
        with open ( os.path.join ( home, 'parameter', 'amber99-residues.ff' ) ) as f:
            for line in f:
                if line.startswith ( 'residue' ):
                    res = line.strip().split()[1]
                    s += line
                    break
                elif line.startswith ( 'chain' ):
                    res = line[-4:-1]

                    if 'nterm' in line:
                        res += 'nt'
                    elif 'cterm' in line:
                        res += 'ct'
                    s += line
                elif line.startswith ( 'bbatom' ):
                    cols = line.strip().split()
                    at = cols[2]
                    s += ' '.join ( cols[:3] )
                    s += ' {0} {0}\n'.format( self.cljs.cljbyres ( res, at ).id )

                elif line.startswith ( 'atom' ):
                    cols = line.strip().split()
                    at = cols[1]
                    s += ' '.join ( cols[:2] )
                    s += ' {0} {0} '.format( self.cljs.cljbyres ( res, conv['backbone'][at] ).id )
                    s += ' '.join ( cols[4:] ) + '\n'

                else:
                    s += line

            for line in f:
                if line.startswith ( 'residue' ):
                    res = line.strip().split()[1]
                    s += line
                elif line.startswith ( 'atom' ):
                    cols = line.strip().split()
                    at = cols[1]
                    s += ' '.join ( cols[:2] )
                    s += ' {0} {0} {1} {1} {2} {2} '.format ( self.cljs.cljbyres ( res, conv[res][at] ).id,
                                                              self.cljs.cljbyres ( res+'nt', conv[res][at] ).id,
                                                              self.cljs.cljbyres ( res+'ct', conv[res][at] ).id )                                                             
                    s += ' '.join ( cols[-3:] ) + '\n'
                else:
                    s += line
        with open ( fname, 'w' ) as f:
            f.write ( s )

class PmsChain:
    def __init__ ( self ):
        self.atoms = {}
        self.rest = []

    def add_atom ( self, name, param_id ):
        self.atoms[name] = param_id

    def __setitem__ ( self, key, val ):
        self.add_atom ( key, val )

    def __getitem__ ( self, key ):
        return self.atoms[key]
            
class PmsResidue:
    def __init__ ( self ):
        self.chains = {}
        self.sidechains = { 'nterm' : {},
                            'cterm' : {},
                            'center' : {} }
        self.rest = []
        
    def add_chain ( self, pos ):
        self.chains[pos] = PmsChain()
        return self.chains[pos]

    def __getitem__ ( self, val ):
        return self.chains[val]

    def add_atom ( self, name, nterm, center, cterm ):
        self.sidechains['nterm'][name] = nterm
        self.sidechains['center'][name] = center
        self.sidechains['cterm'][name] = cterm
            
class ResiduesFile:
    def __init__ ( self, fname = None ):
        self.reses = {}
        if fname != None:
            self.parse ( fname )

    def set_res ( self, name ):
        try:
            return self.reses[name]
        except KeyError: #create new res if one does not already exist
            res = PmsResidue ()
            self.reses[name] = res
            return res
        
    def parse ( self, fname ):
        with open ( fname ) as f:
            for line in f:
                if line.startswith ( 'residue' ):
                    mode = 'res'
                    res_name = line.split()[1].strip()
                    res = self.set_res ( res_name )
                    
                elif line.startswith ( 'chain' ):
                    mode = 'chain'
                  
                    res_name = line[-4:-1]
                    res = self.set_res ( res_name )
                    
                    pos = line[9:-4]
                    chain = res.add_chain ( pos )

                elif line.startswith ( ( 'atom', 'bbatom' ) ):
                    cols = line.split ()
                    if mode == 'chain':
                        chain[cols[2]] = int ( cols[3] )
                    elif mode == 'res':
                        res.add_atom ( cols[1], int ( cols[2] ), int ( cols[4] ), int ( cols[6] ) )
                        

                elif line.startswith ( ( 'zmat', 'bond', 'angle', 'dihedral', 'info', 'backbone' ) ):
                    if mode == 'chain':
                        chain.rest.append ( line )
                    elif mode == 'res':
                        res.rest.append ( line )
