parfile $PROTOMSHOME/parameter/amber14SB.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber14SB-residues.ff
parfile $PROTOMSHOME/parameter/gaff16.ff
parfile l11.tem
solute1 l11.pdb
outfolder out_gas
streamheader off
streamdetail detail
streamenergy energy
streamwarning warning
streaminfo info
streamfatal fatal
streamresults results
streamaccept accept
cutoff 10.0
feather 0.5
temperature 25.0
ranseed 1742449
chunk transrot 1 0.0 0.0
chunk singlepoint
chunk soluteenergy 1
chunk fakesim
chunk results write
