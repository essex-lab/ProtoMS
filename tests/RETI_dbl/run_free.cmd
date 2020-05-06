parfile $PROTOMSHOME/parameter/amber14SB.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber14SB-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
parfile eth-meo.tem
solute1 ethane.pdb
solute2 methanol.pdb
solvent1 ethane_box.pdb
outfolder out_free
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
chunk transrot 2 0.0
printfe mbar
dualtopology1 1 2 synctrans syncrot
softcore1 solute 1
softcore2 solute 2
softcoreparams coul 1 delta 0.2 deltacoul 2.0 power 6 soft66
dlambda 0.001
lambdare 2 0.000 0.330 0.670 1.000
dump 1 results write results
dump 1 results writeinst results_inst
dump 1 pdb all solvent=all file=all.pdb standard
dump 1 restart write restart
dump 1 averages reset
chunk equilibrate 0 solvent=854 protein=0 solute=143 volume=3
chunk simulate 10 solvent=854 protein=0 solute=143 volume=3
