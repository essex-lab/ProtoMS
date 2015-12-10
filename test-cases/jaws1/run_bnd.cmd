parfile $PROTOMSHOME/parameter/amber99.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber99-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
parfile fragment.tem
protein1 protein_scoop.pdb
solute1 fragment.pdb
solvent1 water_clr.pdb
outfolder out
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
#  JAWS-1 specific parameters
jaws1 0
parfile $PROTOMSHOME/data/gcmc_wat.tem
grand1 jaws1_wat.pdb
originx 59.774
originy 31.007
originz 20.931
x 13.02
y 8.188
z 8.402
#  End of JAWS specific parameters
dump 10 results write results
dump 10 pdb all solvent=all file=all.pdb standard
dump 10 restart write restart
dump 10 averages reset
chunk equilibrate 0 solvent=376 protein=124 solute=1 theta=333 gcsolute=167
chunk simulate 100 solvent=376 protein=124 solute=1 theta=333 gcsolute=167