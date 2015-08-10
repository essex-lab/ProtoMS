parfile $PROTOMSHOME/parameter/amber99.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber99-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
parfile dcb.tem
protein1 protein_scoop.pdb
solute1 dcb.pdb
solvent1 water.pdb
outfolder out_bnd
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
chunk equilibrate 100 solvent=841 protein=157 solute=1 volume=0
chunk pdb all solvent=all file=equil_bnd.pdb standard
