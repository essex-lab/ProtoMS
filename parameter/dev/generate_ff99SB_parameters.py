import os, sys
home = os.environ['PROTOMSHOME']
sys.path.append ( os.path.join ( home, 'tools' ) )
import AmbParam, prmtop

ambParm = AmbParam.AmbParameterSet ()

amb_home = os.environ['AMBERHOME']
#must read in .dat first followed by .frcmod then the .in files
ambParm.read_dat ( '%s/dat/leap/parm/parm99.dat' % amb_home )
ambParm.read_frcmod ( '%s/dat/leap/parm/frcmod.ff99SB' % amb_home )

ambParm.read_in ( '%s/dat/leap/prep/oldff/all_amino94.in' % amb_home )
ambParm.read_in ( '%s/dat/leap/prep/oldff/all_aminoct94.in' % amb_home,
                  term = 'cterm' )
ambParm.read_in ( '%s/dat/leap/prep/oldff/all_aminont94.in' % amb_home,
                  term = 'nterm' )



ambParm.fix_diffs ()
s = ambParm.write_protoms_ff( 'amber99SB.ff' )
s = ambParm.write_protoms_residues( 'amber99SB-residues.ff' )
