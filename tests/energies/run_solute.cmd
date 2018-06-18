parfile $PROTOMSHOME/parameter/amber14SB.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber14SB-residues.ff
parfile brd.tem
parfile t3p.tem
solute1 brd.pdb
solute2 cld.pdb
solute3 soft_solute_t3p.pdb
outfolder out_solute
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

dualtopology1 1 2 synctrans syncrot
lambdas 0.5 0.75 0.25
chunk fakesim
chunk results write results
