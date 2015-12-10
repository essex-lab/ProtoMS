parfile $PROTOMSHOME/parameter/amber99.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber99-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
parfile ethtmeo_ele.tem
solute1 ethane.pdb
solvent1 ethane_box.pdb
outfolder out_ele_free
streamheader off
streamdetail off
streamwarning warning
streaminfo info
streamfatal fatal
streamresults results
streamaccept accept
cutoff 10.0
feather 0.5
temperature 25.0
ranseed 100000
boundary solvent
pressure 1
chunk transrot 1 0.0
printfe mbar
dlambda 0.001
lambdare 2 0.000 0.330 0.670 1.000
dump 1 results write results
dump 1 results writeinst results_inst
dump 1 pdb all solvent=all file=all.pdb standard
dump 1 restart write restart
dump 1 averages reset
chunk equilibrate 0 solvent=981 protein=0 solute=16 volume=3
chunk simulate 10 solvent=981 protein=0 solute=16 volume=3