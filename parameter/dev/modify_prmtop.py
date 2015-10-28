import os, sys, argparse
home = os.environ['PROTOMSHOME']
sys.path.append ( os.path.join ( home, 'tools' ) )
import prmtop

parser = argparse.ArgumentParser ( description = "Script to modify amber topology files to omit dihedral energetic terms not calculated by ProtoMS and rescale atomic charges. The resultant topology file can be used in validation of parameter sets used in ProtoMS."  )
parser.add_argument ( '-i', '--infile', required = True, help = 'Name of input topology file to be modified.' )
parser.add_argument ( '-o', '--outfile', required = True, help = 'Name of output topology file.' )
args = parser.parse_args()

    
#parse input topology file and pull out those sections we'll need
prm = prmtop.prmtop ( args.infile )
s = prm.sec_dict['DIHEDRALS_INC_HYDROGEN']
s2 = prm.sec_dict['DIHEDRALS_WITHOUT_HYDROGEN']
names = prm.sec_dict['ATOM_NAME'].vals
res_point = prm.sec_dict['RESIDUE_POINTER']
res_names = prm.sec_dict['RESIDUE_LABEL']
charges = prm.sec_dict['CHARGE']

#define dihedrals to be removed
sidechain_impropers = { 'ARG' : ( ('NH2', 'CZ', 'NH1', 'NE'),
                                  ('CD', 'CZ', 'NE', 'HE'),
                                  ('CZ', 'HH11', 'NH1', 'HH12'),
                                  ('CZ', 'HH21', 'NH2', 'HH22') ),
                        'ASN' : ( ('OD1', 'CG', 'ND2', 'CB'),
                                  ('CG', 'HD21', 'ND2', 'HD22'),
                                  ('CG', 'HD21', 'ND2', 'HD22') ),
                        'ASP' : ( ('OD2', 'CG', 'OD1', 'CB'), ),
                        'GLN' : ( ('OE1', 'CD', 'NE2', 'CG'),
                                  ('CD', 'HE21', 'NE2', 'HE22') ),
                        'GLU' : ( ('OE2','CD','OE1','CG'), ),
                        'HIS' : ( ('CB', 'CD2', 'CG', 'ND1'),
                                  ('HD2', 'CD2', 'NE2', 'CG'),
                                  ('CE1', 'CD2', 'NE2', 'HE2'),
                                  ('HE1', 'CE1', 'NE2', 'ND1'),
                                  ('CG', 'CE1', 'ND1', 'HD1') ) ,
                        'PHE' : ( ('CD1', 'CD2', 'CG', 'CB'),
                                  ('HD2', 'CD2', 'CE2', 'CG'),
                                  ('CZ', 'CD2', 'CE2', 'HE2'),
                                  ('CE1', 'CE2', 'CZ', 'HZ'),
                                  ('HE1', 'CE1', 'CZ', 'CD1'),
                                  ('HD1', 'CD1', 'CE1', 'CG') ), 
                        'TRP' : ( ('CD2', 'CB', 'CG', 'CD1'),
                                  ('CZ3', 'CD2', 'CE3', 'HE3'),
                                  ('CH2', 'CE3', 'CZ3', 'HZ3'),
                                  ('HH2', 'CH2', 'CZ3', 'CZ2'),
                                  ('CE2', 'CH2', 'CZ2', 'HZ2'),
                                  ('CD1', 'CE2', 'NE1', 'HE1'),
                                  ('HD1', 'CD1', 'NE1', 'CG') ),
                        'TYR' : ( ('CE1', 'CE2', 'CZ', 'OH'),
                                  ('CD1', 'CD2', 'CG', 'CB'),
                                  ('HD2', 'CD2', 'CE2', 'CG'),
                                  ('CZ', 'CD2', 'CE2', 'HE2'),
                                  ('HE1', 'CE1', 'CZ', 'CD1'),
                                  ('HD1', 'CD1', 'CE1', 'CG') ) }

sidechain_rings = { 'HIS' : ( ('CG', 'CD2', 'NE2', 'CE1'),
                              ('NE2', 'CE1', 'ND1', 'CG'),
                              ('ND1', 'CE1', 'NE2', 'CD2'),
                              ('NE2', 'CD2', 'CG', 'ND1'),
                              ('CD2', 'CG', 'ND1', 'CE1') ),
                    'PHE' : ( ('CG', 'CD2', 'CE2', 'CZ'),
                              ('CD1', 'CE1', 'CZ', 'CE2'),
                              ('CD2', 'CE2', 'CZ', 'CE1') ),
                    'PRO' : ( ('N', 'CD', 'CG', 'CB'),
                              ('CB', 'CA', 'N', 'CD'),
                              ('CG', 'CD', 'N', 'CA') ),
                    'TRP' : ( ('CD2', 'CE3', 'CZ3', 'CH2'),
                              ('CZ2', 'CH2', 'CZ3', 'CE3'),
                              ('CE2', 'CD2', 'CE3', 'CZ3'),
                              ('NE1', 'CD1', 'CG', 'CD2'),
                              ('CE2', 'CD2', 'CG', 'CD1'),
                              ('CD2', 'CE2', 'NE1', 'CD1'),
                              ('CG', 'CD1', 'NE1', 'CE2'),
                              ('CG', 'CD2', 'CE2', 'NE1') ) }
sidechain_rings['TYR'] = sidechain_rings['PHE']

backbone_impropers = ( ( 'O', 'C', 'N', 'CA' ),
                       ( 'C', 'CA', 'N', 'H' ),
                       ( 'C', 'CD', 'N', 'CA' ) )


def get_res_info ( at ):
    """Given an at id find the corresponding residue number 
    and name information from the prmtop file"""
    for i in xrange ( len ( res_point.vals ) -1, -1, -1 ):
        if at >= res_point.vals[i]:
            return i + 1, res_names.vals[i].strip()

#build list of dihedrals to remove
to_pop = []
for i, chnk in enumerate ( prmtop.chunks ( s.vals + s2.vals, 5 ) ):
    ids = map ( prmtop.ptoi,
                map ( abs, chnk[:4] ) )
    dih_ats = tuple ( names[ j - 1 ].strip() for j in ids )
    reses = [ get_res_info ( j )[1] for j in ids ]
    res_ids = [ get_res_info ( j )[0] for j in ids ]    
    
    #flipping rules to make atom order consistent between amb and pms
    if dih_ats[3] < dih_ats[0]:
        dih_ats = dih_ats[::-1]
        reses = reses[::-1]
        res_ids = res_ids[::-1]
    if dih_ats[2] < dih_ats[1]:
        dih_ats = dih_ats[::-1]
        reses = reses[::-1]
        res_ids = res_ids[::-1]
    
    #Treat all HIS based residues the same for simplicity
    reses = [ 'HIS' if res in ( 'HIE', 'HID', 'HIP' ) else res for res in reses ]        
    
    try:
        #As we are interested in sidechains we can safely assume that the
        #residue name of the first atom is the same for all
        for j in sidechain_impropers[reses[0]]: 
            if j == dih_ats and chnk[2] < 0:
                to_pop.append ( i )
    except KeyError:
        pass

    try:
        for j in sidechain_rings[reses[0]]:
            if j == dih_ats:
                to_pop.append ( i )
    except KeyError:
        pass
            
    #C terminus improper
    if dih_ats == ('OXT','C','O','CA'): 
        to_pop.append ( i  )

    if dih_ats in backbone_impropers and chnk[2] < 0:
        to_pop.append ( i )

#Now remove dihedrals from prmtop file
to_pop = sorted ( set ( to_pop ) )
for i in to_pop[::-1]:
    try:
        for k in xrange ( 5 ):
            s.vals.pop ( 5*i )
    except IndexError:
        for k in xrange ( 5):
            s2.vals.pop ( 5*i - len ( s.vals ) )

#Adjust values in pointer sections to be consistent
ps = prm.sec_dict [ 'POINTERS' ]
ps.vals[6] = len ( s.vals ) / 5
ps.vals[7] = len ( s2.vals ) / 5
ps.vals[14] = ps.vals[7]

#Now adjust charges
#The charges in prmtop files are in internal amber units
#These are produced by scaling atomic charges by 18.2223
#The truncation of this factor at 4 d.p. produces
#discrepancies in electrostatic energies between ProtoMS
#and AMBER when compounded over a large number of interactions.

#To reduce this difference we here readjust the charges such that
#the conversion factor is the more correct value of 18.2222615
charges.vals = [ i / 18.2223 * 18.222615 for i in charges.vals ]

with open ( args.outfile, 'w' ) as f:
    f.write ( str ( prm ) )



    
####Old testing stuff below here#####




# #get list of remaining dihedrals after popping
# remain_list = []
# for i, chnk in enumerate ( prmtop.chunks ( s.vals + s2.vals, 5 ) ):
#     ids = map ( prmtop.ptoi,
#                 map ( abs,
#                     map ( int, chnk[:4] ) ) )
#     dih_ats = tuple ( names[ j - 1 ].strip() for j in ids )
#     reses = [ get_res_info ( j )[1] for j in ids ]
#     res_ids = [ get_res_info ( j )[0] for j in ids ]    
    
#     #flipping rules to make atom order consistent between amb and pms
#     if dih_ats[3] < dih_ats[0]:
#         dih_ats = dih_ats[::-1]
#         reses = reses[::-1]
#         res_ids = res_ids[::-1]
#     if dih_ats[2] < dih_ats[1]:
#         dih_ats = dih_ats[::-1]
#         reses = reses[::-1]
#         res_ids = res_ids[::-1]
#     remain_list.append ( dih_ats + tuple ( res_ids ) )



# # for i in to_pop:
# #     try:
# #         remain_list.index ( dih_list[i] )
# #     except ValueError:
# #         continue
# #     print i

# for i in remain_list:
#     if i not in dih_list:
#         print i


# pop_list = [ dih_list[i] for i in to_pop ]

# t = ( 'C','CA','CB','CG', 1, 1, 1, 1 )
# t2 = dih_list[to_pop[0]]
# print t in dih_list
# print t in pop_list
# print t in remain_list
# print
# print t2 in dih_list
# print t2 in pop_list
# print t2 in remain_list
# print

# for i in pop_list:
#     if i in remain_list:
#         print i


# for i in remain_list + pop_list:
#     if i not in dih_list:
#         print i
# for i in dih_list:
#     if i not in remain_list + pop_list:
#         print i

# for i in remain_list:
#     if remain_list.count ( i ) > dih_list.count ( i ):
#         print i

        
