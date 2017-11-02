parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/gaff16.ff
parfile t3p.ff
parfile brd.tem
protein1 soft_protein.pdb
solute1 soft_solute.pdb
solute2 soft_solute.pdb
outfolder out_protein
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

