*****************
protoms.py
*****************

This program is used to setup a ProtoMS simulation. It was made with usability at highest priority. The only input that should be necessary is a couple of prepared PDB files containing the molecules one would like to simulate. 

The program will create force field for small molecules, setup the protein and solvate the prepared system. At the moment it can setup the following types of simulations:

* Equilibration
* Sampling
* Dual-topology free energy
* Single-topology free energy
* Grand Canonical Monte Carlo (GCMC)
* Just Add Waters, stage 1 and 2 (JAWS-1, JAWS-2)

The program will create files and inputs based on experience that should work in most situations. However there might be situations where the created settings are not appropriate. One can then use individual tools to make a more custom setup, see this . One might also have to edit the files manually.


**Syntax:**

  ``protoms.py [-s none|equilibration|sampling|dualtopology|singletopology|gcmc|jaws1|jaws2] [-f folder1 folder2] [-p protein.pdb] [-o scoop.pdb] [-l lig1.pdb lig2.pdb ...] [-t template1 template2 ...] [-w water.pdb] [-c cmdfile] [-r nrepeats | prefix] [--outfolder folder] [--atomnames namefile] [--watmodel tip4p|tip3p] [--waterbox watbox] [--charge charge1 charge2] [--singlemap mapfile]  [--center cent] [--innercut icut] [--outercut ocut] [--flexin sidechain|flexible|rigid] [--flexout sidechain|flexible|rigid] [--scooplimit N] [--capradius radius] [--lambdas nlambdas | lambda1 lambda2 ...] [--adams B1 B2 ...] [--jawsbias bias] [--gcmcwater wat.pdb] [--gcmcbox box.pdb] [--nequil N] [--nprod N] [--dumpfreq N] [--absolute] [--dovacuum] [--testrun]``


*  ``-s none|equilibration|sampling|dualtopology|singletopology|gcmc|jaws1|jaws2`` = the type of simulation to perform
    optional, default = ``none``
*  ``-f folder1 folder2`` = name of folders to search for input files
    optional, no default
* ``-p protein.pdb`` = the name of the protein PDB file
    optional, no default
* ``-o scoop.pdb`` = the name of a protein scoop PDB file
    optional, no default    
* ``-l lig1.pdb lig2.pdb ...`` = the name(s) of PDB file(s) containing ligand(s)
    optional, no default
* ``-t template1 template2 ...`` = the name(s) of ProtoMS template file(s) that needs to be loaded
    optional, no default
* ``-w water.pdb`` = the name of a PDB file with bulk water for the protein
    optional, no default
* ``-c cmdfile`` = the prefix for the created ProtoMS command file
    optional, default = ``run``
* ``-r nrepeats | prefix`` = setup independent repeats of the simulation
    optional, default = 1
    ``nrepeat`` = repeats a created from 1 to ``nrepeat``
    ``prefix`` = a single repeat is created, but prefix is appended to folders and files
* ``--outfolder folder`` = the ProtoMS output folder
    optional, default = ``""`` (empty string)
* ``--atomnames namefile`` = the name of file containing conversion instructions
    optional, no default
    if not given, takes the one in ``$PROTOMSHOME/data``
* ``--watermodel tip4p|tip3p`` = the water model to use
    optional, default = ``tip4p``
* ``--waterbox watbox`` = the name a of a PDB file with a pre-equilibrated water box
    optional, no default
    if not given, takes one in ``$PROTOMSHOME/data``
* ``--charge charge1 charge2`` ... = the charges of the ligands
    optional, default = 0
* ``--singlemap mapfile`` = the correspondence map for single-topology setup
    optional, no default
* ``--center cent`` = the centre of the scoop
    optional, default = 0.0,0.0,0.0
* ``--innercut icut`` == the inner region cut-off in Angstroms
    optional, default = 16.9 A
* ``--outercut ocut`` == the outer region cut-off in Angstroms
    optional, default = 20.0 A
* ``--flexin sidechain|flexible|rigid`` = determine the flexibility of the inner region
    optional, default = ``flexible``
    ``sidechain`` = only the sidechains will be sampled in the simulation
    ``flexible`` = both sidechain and backbone will be sampled in the simulation
    ``rigid`` = no residues will be sampled
* ``--flexout sidechain|flexible|rigid`` = determine the flexibility of the outer region
    optional, default = ``sidechain``
    ``sidechain`` = only the sidechains will be sampled in the simulation
    ``flexible`` = both sidechain and backbone will be sampled in the simulation
    ``rigid`` = no residues will be sampled
* ``--scooplimit N`` = the minimum removed number of residues in a scoop
    optional, default = 10
* ``--capradius radius`` = the radius of the droplet solvating the protein
    optional, default = 30  
* ``--lambdas nlambdas | lambda1 lambada2`` ... = specification of :math:`\lambda`; space for free energy calculations
    optional, default = 16
    if a single value is given, this number of :math:`\lambda`-values is created uniformly from 0 to 1
    if a list of values are given, this is the :math:`\lambda`-values to use
* ``--adams B1 B2`` ... = the Adams parameter for GCMC
    optional, default = 0
* ``--jawsbias bias`` = the bias to apply in JAWS-2 simulations
    optional, default = 0
* ``--gcmcwater wat.pdb`` = the name of a PDB file with reservoir waters for GCMC and JAWS-1
    optional, no default
* ``--gcmcbox box.pdb`` = the name of a PDB file with GCMC or JAWS-1 simulation box dimension
    optional, no default
* ``--nequil N`` = the number of equilibration moves
    optional, default = 5E6
* ``--nprod N`` = the number of production moves
    optional, default = 40E6
* ``--dumpfreq N`` = the frequency with which output is written to disc
    optional, default = 1E5
* ``--absolute`` = turns *on* the setup of absolute free energies
    optional, default = off          
* ``--dovacuum`` = turns *on* the setup of vacuum simulation
    optional, default = off
* ``--testrun`` = turns *on* the setup of a short simulations appropiate for tests
    optional, default = off


**Examples:**

:: 

  protoms.py
  protoms.py.py -s sampling -l lig1.pdb --dovacuum --testrun
  protoms.py -s dualtopology -l lig1.pdb lig2.pdb -p protein.pdb
  protoms.py -s dualtopology -l lig1.pdb --absolute
  protoms.py -s gcmc -p protein.pdb --adams -4 -2 0 2 4 6

**Notes:**

The program will try to locate previously created files for the protein and ligand in the current working directory or any folder specified with the ``-f`` flag. For ligands the program will replace ``.pdb`` with the appropriate ending, such as ``.prepi`` for Amber prepi files and ``.tem`` for ProtoMS template files.

Starting with just the PDB-files of the ligand(s) and the protein, the program will create the following files in the same folder as those PDB-files

* ``lig.prepi`` = the z-matrix and atom types of the ligand in Amber format
* ``lig.frcmod`` = additional parameters not in GAFF
* ``lig.zmat`` = the z-matrix of the ligand used to sample it in the MC simulation
* ``lig.tem`` = the complete template (force field) file for the ligand in ProtoMS format
* ``li1-li2.tem`` = the combined template file of all ligands 
   the filename is a combination of the residue name of all ligands   
* ``lig_box.pdb`` = the box of water solvating the ligand
* ``protein_scoop.pdb`` = the truncated protein structure
* ``protein_pms.pdb`` = the original protein structure with ProtoMS naming convention
   if the scoop removes to few residues, this file be created instead
* ``water.pdb`` = the cap of water solvating the protein system

In addition, for dual-topology simulations the following files are created: :

* ``lig1_dummy.pdb`` = the dummy particle that the ligand will be perturbed to
   only created if the --absolute flag is set 

In addition, for single-topology simulations the following files are created:

* ``li1-li2_ele.tem`` = the ProtoMS template file for electrostatic single-topology perturbation 
* ``li1-li2_vdw.tem`` = the ProtoMS template file for van der Waals single-topology perturbation       
* ``settings.singlemap`` = the created correspondance map for single topology  
    only named like this if the --singlemap argument is not set  


In addition, for GCMC / JAWS-1 simulations the following files are created:


* ``gcmc_box.pdb`` / ``jaws1_box.pdb`` = the GCMC / JAWS-1 simulation box 
* ``gcmc_wat.pdb`` = the GCMC / JAWS-1 reservoire waters       
* ``water_clr.pdb`` = the cap of water solvating the protein system, cleared from the GCMC / JAWS-1 simulation box


In addition, for JAWS-2 simulations the following files are created:

* ``jaws2_watN.pdb`` = the JAWS-2 water
    each of the water given with the ``--gcmc_water`` flag will be written to an individual file 
* ``jaws2_notN.pdb`` = the rest of the JAWS-2 water
* ``water_clr.pdb`` = the cap of water solvating the protein system, cleared from the GCMC / JAWS-1 simulation box


It will create at most three ProtoMS command files, one for the protein simulation, one for the ligand simulation and one for the gas-phase simulation. These can be used to run ProtoMS, e.g. ::

  $PROTOMS/protoms3 run_free.cmd

**Prerequisites:**

The program assumes that both the ligand and the protein is prepared before. This includes for instance protonation. One can read more about setting up ligands and proteins here.

The progam requires AmberTools to make force field for small molecules.

