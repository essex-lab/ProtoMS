import os
import sys

home = os.environ["PROTOMSHOME"]
sys.path.append(os.path.join(home, "tools"))
import AmbParam

ambParm = AmbParam.AmbParameterSet(1001)
# ambParm2 = AmbParam.AmbParameterSet ()

amb_home = os.environ["AMBERHOME"]
# Gaff ff's only have a .dat file associated
# When development of ff is complete, also possible to run this for Gaff v2 below
ambParm.read_dat("%s/dat/leap/parm/gaff.dat" % amb_home)
# ambParm2.read_dat ( '%s/dat/leap/parm/gaff2.dat' % amb_home )


# fix_gaff_diffs removes the hw-ow bond defined for the Amber fast water routine
ambParm.fix_gaff_diffs()
s = ambParm.write_protoms_types("gaff16.types")
s = ambParm.write_protoms_ff("gaff16.ff")
# s = ambParm2.write_protoms_types( 'gaff-v2.types' )
# s = ambParm2.write_protoms_ff( 'gaff-v2.ff' )
