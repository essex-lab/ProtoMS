parfile $PROTOMSHOME/parameter/amber14SB.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber14SB-residues.ff
parfile brd.tem
solute1 soft_solute.pdb
solute2 soft_solute.pdb
outfolder out_gcsolute
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
pdbparams on
gcmc 1
parfile $PROTOMSHOME/data/gcmc_tip3p.tem
grand1 soft_gcsolute.pdb
potential 0

dualtopology1 1 2 synctrans syncrot
lambdas 0.5 0.75 0.25
chunk fakesim
chunk results write results
