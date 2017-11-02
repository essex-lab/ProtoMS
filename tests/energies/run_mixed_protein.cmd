parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/gaff16.ff
parfile t3p.ff
parfile brd.tem
protein1 soft_protein.pdb
solute1 soft_solute.pdb
solute2 soft_solute.pdb
outfolder out_mixed_protein
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
softcore1 solute 1 atom 1
softcore2 solute 2 atom 1
softcoreparams coul 1 delta 0.2 deltacoul 2.0 power 0 soft66
lambdas 0.5 0.75 0.25
chunk fakesim
chunk results write results

