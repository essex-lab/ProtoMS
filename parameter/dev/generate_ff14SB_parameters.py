import os
import sys

home = os.environ["PROTOMSHOME"]
sys.path.append(os.path.join(home, "tools"))
import AmbParam

ambParm = AmbParam.AmbParameterSet()

amb_home = os.environ["AMBERHOME"]
# must read in .dat first followed by .frcmod then the .in files
ambParm.read_dat("%s/dat/leap/parm/parm10.dat" % amb_home)
ambParm.read_frcmod("%s/dat/leap/parm/frcmod.ff14SB" % amb_home)

ambParm.read_in("%s/dat/leap/prep/amino12.in" % amb_home)
ambParm.read_in("%s/dat/leap/prep/aminoct12.in" % amb_home, term="cterm")
ambParm.read_in("%s/dat/leap/prep/aminont12.in" % amb_home, term="nterm")


ambParm.fix_diffs()
s = ambParm.write_protoms_ff("amber14SB.ff")
s = ambParm.write_protoms_residues("amber14SB-residues.ff")
