import os, sys
home = os.environ['PROTOMSHOME']
sys.path.append ( os.path.join ( home, 'tools' ) )
from convertatomnames import read_convfile
import simulationobjects as sim
import prmtop
conv = read_convfile ( os.path.join ( home, 'data', 'atomnamesmap2.dat' ),
                       inmode = 'protoms', outmode = 'amber' )

pms_parms = sim.ParameterSet ( 'dihedral' )
pms_parms.read ( 'tmp.ff' )

pms_reses = sim.ResiduesFile ( 'tmp-residues.ff' )

prm = prmtop.prmtop ( 'test/prmtop_test' )


s = prm.sec_dict['DIHEDRALS_INC_HYDROGEN']
s2 = prm.sec_dict['DIHEDRALS_WITHOUT_HYDROGEN']
names = prm.sec_dict['ATOM_NAME'].vals
amb_types = prm.sec_dict['AMBER_ATOM_TYPE'].vals
ks = prm.sec_dict['DIHEDRAL_FORCE_CONSTANT']
phase = prm.sec_dict['DIHEDRAL_PHASE']

amb_Parm = sim.AmbParameterFile ()

#must read in .dat first followed by .in files in any order
amb_Parm.read_dat ( '/home/chris/Software/amber14/dat/leap/parm/parm99.dat' )
amb_Parm.read_in ( '/home/chris/Software/amber14/dat/leap/prep/oldff/all_amino94.in' )
amb_Parm.read_in ( '/home/chris/Software/amber14/dat/leap/prep/oldff/all_aminoct94.in',
                  term = 'cterm' )
amb_Parm.read_in ( '/home/chris/Software/amber14/dat/leap/prep/oldff/all_aminont94.in',
                  term = 'nterm' )
amb_Parm.fix_diffs()


#build dictionary of amber dihedral counts
amb_dihs = {}
for i in prmtop.chunks ( s.vals + s2.vals, 5 ):
    # if dihedral prefactor is zero then don't count it
    if float ( ks.vals[ int ( i[4]  ) - 1 ] ) == 0.0:
        continue
    
    ids = map ( prmtop.ptoi,
                map ( abs,
                    map ( int, i[:4] ) ) )
    dih_ats = tuple ( names[ j - 1 ].strip() for j in ids )

    #flipping rules to make atom order consistent between amb and pms
    if dih_ats[3] < dih_ats[0]:
        dih_ats = dih_ats[::-1]
    if dih_ats[2] < dih_ats[1]:
        dih_ats = dih_ats[::-1]        

    try:
        amb_dihs[ dih_ats ] += 1
    except KeyError:
        amb_dihs[ dih_ats ] = 1

    # if dih_ats == ('O', 'C', 'N', 'CA'):
    #     # if dih_ats == ('O', 'C', 'N', 'CA'):
    #     print ks.vals[ int ( i[4] )  - 1 ], dih_ats, ids, [ amb_types[ prmtop.ptoi ( abs ( int ( j ) ) ) - 1] for j in i[:4]  ], i

pms_dihs = {}
# with open ( 'test2/out_bnd/detail' ) as f:
with open ( 'test/tmp_out' ) as f:
    for line in f:
        if line.startswith ( ' TAG' ):
            cols = line.upper().strip().split()
            ats = cols[1:17:4]
            reses = cols[2:19:4]
            res_nos = [ i.strip(')') for i in cols[4:21:4] ]

            # print cols
            # break
            
            amb_ats = []
            for i, j in zip ( ats, reses ):
                try:
                    amb_ats += [ conv[j.upper()][i.upper()] ]
                except KeyError:
                    amb_ats += [ conv['backbone'][i.upper()] ]

            #fix res names so that we are looking at correct templates for terminals
            reses = [ i + 'nt' if j == '1' else i for i, j in zip ( reses, res_nos )]
            reses = [ i + 'ct' if j == '157' else i for i, j in zip ( reses, res_nos ) ]

            #flip rules as above
            if amb_ats[3] < amb_ats[0]:
                amb_ats = amb_ats[::-1]
                reses = reses[::-1]
            if amb_ats[2] < amb_ats[1]:
                amb_ats = amb_ats[::-1]
                reses = reses[::-1]

            #get atom types and look up parameters
            at_types = [ amb_Parm.cljs.cljbyres (  i, j ).type for i, j in zip ( reses, amb_ats ) ]
            dih =  pms_parms.get_params ( at_types )
            nterms = len (  [ i for i in dih.terms if i.k1 != 0.0 ] ) # get non-zero dihedral terms

            # if amb_ats == [ 'C', 'CA', 'N', 'H' ]:
            #     break
            
            # if amb_ats == ['O', 'C', 'N', 'CA']:
            #     # print at_types, dih.terms, nterms
            #     print at_types, dih.terms
            #     # break
            # for t in dih.terms:
            #     print t.k1, t.k2, t.k3, t.k4
            # print
            # if nterms == 0:
            #     print at_types, dih.ats
            try:
                pms_dihs[ tuple ( amb_ats ) ] +=   nterms
            except KeyError:
                pms_dihs[ tuple ( amb_ats ) ] = nterms
                

print(sum(pms_dihs.values()), sum(amb_dihs.values()))


all_ats = set ( pms_dihs.keys() + amb_dihs.keys() )
                
for i in all_ats:
    try:
        ambn = amb_dihs[i]
    except KeyError:
        ambn = 0
    try:
        pmsn = pms_dihs[i]
    except KeyError:
        pmsn = 0

    if ambn != pmsn:
        print i, pmsn, ambn
                
