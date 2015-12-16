parfile $PROTOMSHOME/parameter/amber99.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber99-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
parfile fragment.tem
protein1 protein_scoop.pdb
solute1 fragment.pdb
solvent1 water.pdb
outfolder out_jaws2-w2
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
#  JAWS-2 specific parameters
jaws2 1
multijaws2 8.000 10.000 12.000 14.000
parfile $PROTOMSHOME/data/gcmc_wat.tem
solvent2 jaws2_not2.pdb
grand1 jaws2_wat2.pdb
originx 59.271
originy 33.21
originz 21.691
x 3.0
y 3.0
z 3.0
#  End of JAWS specific parameters
dump 10 results write results
dump 10 pdb all solvent=all file=all.pdb standard
dump 10 restart write restart
dump 10 averages reset
chunk equilibrate 0 solvent=376 protein=124 solute=1 sample=333 gcsolute=167
chunk simulate 100 solvent=376 protein=124 solute=1 sample=333 gcsolute=167
