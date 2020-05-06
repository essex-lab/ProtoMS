parfile $PROTOMSHOME/parameter/dev/tmp.ff
parfile $PROTOMSHOME/parameter/dev/tmp-residues.ff
protein1 protein_pms.pdb
outfolder out_bnd
streamheader off
streamdetail detail
streamwarning warning
streaminfo info
streamfatal fatal
streamresults results
streamaccept accept
cutoff 500.0
feather 0.0
temperature 25.0
ranseed 7821995
boundary solvent
pdbparams on

chunk fixresidues 1 none
chunk fixbackbone 1 none
chunk fakesim
chunk results write results
chunk pdb all solvent=all file=test.pdb standard
