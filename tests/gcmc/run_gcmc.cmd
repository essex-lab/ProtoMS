parfile $PROTOMSHOME/parameter/amber14SB.ff
parfile $PROTOMSHOME/parameter/solvents.ff
parfile $PROTOMSHOME/parameter/amber14SB-residues.ff
parfile $PROTOMSHOME/parameter/gaff14.ff
protein1 protein.pdb
solvent1 water.pdb
outfolder out_gcmc
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
#  GCMC specific parameters
gcmc 0
parfile $PROTOMSHOME/data/gcmc_tip4p.tem
grand1 gcmc_wat.pdb
multigcmc 200 19.000 20.000 
originx 30.024
originy 1.952
originz 8.033
x 5.3
y 4.731
z 4.645
#  End of GCMC specific parameters
dump 100 results write results
dump 100 pdb all solvent=all file=all.pdb standard
dump 100 restart write restart
dump 100 averages reset
chunk equilibrate 0 solvent=0 protein=0 solute=0 insertion=333 deletion=333 gcsolute=333
chunk equilibrate 0 solvent=440 protein=60 solute=0 insertion=167 deletion=167 gcsolute=167
chunk simulate 1000 solvent=440 protein=60 solute=0 insertion=167 deletion=167 gcsolute=167
