parfile $PROTOMSHOME/parameter/amber99.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber99-residues.ff
parfile $PROTOMSHOME/parameter/gaff.ff
parfile ethtmeo_comb.tem
solute1 ethane.pdb
outfolder out_single
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
ranseed 669165
chunk transrot 1 0.0 0.0
printfe mbar
lambda 0.0 0.0 0.0
dump 1 pdb all solvent=all file=all.pdb standard
chunk simulate 1 solvent=0 protein=0 solute=1000 volume=0
chunk lambda 0.25
chunk simulate 1 solvent=0 protein=0 solute=1000 volume=0
chunk lambda 0.50
chunk simulate 1 solvent=0 protein=0 solute=1000 volume=0
chunk lambda 0.75
chunk simulate 1 solvent=0 protein=0 solute=1000 volume=0
chunk lambda 1.0
chunk simulate 1 solvent=0 protein=0 solute=1000 volume=0
