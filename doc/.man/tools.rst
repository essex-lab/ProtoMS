*************
Tools
*************

In the ``$PROTOMSHOME/tools`` folder we have collect a range of useful scripts to setup and analyse ProtoMS simulations. Many of them are used by the ``protoms.py`` setup script. In this page we have collected the documentation for these tools with the user as a focus. Developers might be interested in looking at the Python code manual in the ``.doc`` folder.



	  
	  
----------------------------
ambertools.py
----------------------------
.. argparse::
   :module: ambertools
   :func: get_arg_parser
   :prog: ambertools.py

**Examples:**

:: 

  ambertools.py -f benzene.pdb
  ambertools.py -f benzene.pdb -n BNZ
  ambertools.py -f benzenamide.pdb -c 0
  ambertools.py -f benzene.pdb toluene.pdb -n BNZ TOL


**Description:**

This tool encapsulate the program ``antechamber`` and ``parmchk`` from the AmberTools suite of programs.

It will produce an Amber prepi-file, containing the z-matrix and atom types of the given solutes, parametrized with the general Amber force field and AM1-BCC charges. It will also produce an Amber frcmod-file with additional parameters not found in the GAFF definition. These files be named after the input ``pdbfile``, replacing the extension ``.pdb`` with ``.prepi`` and ``.frcmod``

The ``antechamber`` and ``parmchk`` program should exist in the system path or the AMBERHOME environment variable should be set correctly.


-----------------------
build_template.py
-----------------------
.. argparse::
   :module: build_template
   :func: get_arg_parser
   :prog: build_template.py

**Examples:**

::

  build_template.py -p benzene.prepi
  build_template.py -p benzene.prepi -f benzene.frcmod 
  build_template.py -p benzene.prepi -f benzene.frcmod -o benzene.template -n BNZ
  build_template.py -p benzene.prepi -f benzene.frcmod -t 1.0 -r 10
 
**Description:**

This tool builds a ProtoMS template file for a solute given an Amber prepi file.

If the solute needs parameters not in the specified GAFF release, they should be supplied with the ``frcmodfile``. 

The tool will automatically make an appropriate z-matrix for Monte Carlo sampling. This works in most situations. However, if something is not working properly with the generated z-matrix, one can be supplied in the ``zmatfile``

The default translational and rotational displacements are based on experience and should be appropriate in most situations.

-----------------------
calc_clusters.py
-----------------------
.. argparse::
   :module: calc_clusters
   :func: get_arg_parser
   :prog: calc_clusters.py

**Examples:**

::

  calc_clusters.py -i all.pdb
  calc_clusters.py -i all.pdb all2.pdb
  calc_clusters.py -i all.pdb -o all_clusters.pdb
  calc_clusters.py -i all.pdb -t complete

**Description:**

This tool cluster molecules from a simulation

It will extract the coordinates of all atoms with name equal to ``atom`` in residues with name equal to ``molecule`` in all input files and cluster them using the selected algorithm.  If no atom is specified, the entire molecule will be clustered. By default this atom and residue name is set to match GCMC / JAWS output with the standard water template.

-----------------------
calc_density.py
-----------------------
.. argparse::
   :module: calc_density
   :func: get_arg_parser
   :prog: calc_density.py

**Examples:**

::

  calc_density.py -i all.pdb
  calc_density.py -i all.pdb all2.pdb
  calc_density.py -i all.pdb -o gcmc_density.dx
  calc_density.py -i all.pdb -r t4p -n o00
  calc_density.py -i all.pdb -p 1.0 -s 1.0
  calc_density.py -i all.pdb -e 0.5 -t gaussian
  calc_density.py -i all.pdb -n 100

**Description:**

This tool discretises atoms on a grid, thereby representing a simulation output as a density. 

It will extract the coordinates of all atoms with name equal to ``atom`` in residues with name equal to ``residue`` in all input files and discretise them on a grid. By default this atom and residue name is set to match GCMC / JAWS output with the standard water template.

The produced density can be visualized with most programs, e.g. ::

  vmd -m all.pdb grid.dx


-----------------------
calc_dg.py
-----------------------
.. argparse::
   :module: calc_dg
   :func: get_arg_parser
   :prog: calc_dg.py

**Examples:**

::

  calc_dg.py -d out_free/
  calc_dg.py -d out_free1/ out_free2/ out_free3/ -l 0.1
  calc_dg.py -d out_free1/ out_free2/ out_free3/ -u 0.9
  calc_dg.py -d out_free1/ out_free2/ out_free3/ -e ti bar
  calc_dg.py -d out_free1/ out_free2/ out_free3/ -e gcap
  calc_dg.py -d out_free1/ out_free2/ out_free3/ --subdir b_-9.700 -e ti bar

**Description:**

This tool calculates free energies using the method of thermodynamic integration (TI), Bennett's Acceptance Ratio (BAR), multi state BAR (MBAR) and grand canonical alchemical perturbation (GCAP).

The program expects that in the ``directory``, ``directory2`` etc. there exists an output folder for each :math:`\lambda`-value, eg. ``lam-0.000`` and ``lam-1.000`` (unless the ``--subdir`` argument is used.)

-----------------------
calc_dg_cycle.py
-----------------------
.. argparse::
   :module: calc_dg_cycle
   :func: get_arg_parser
   :prog: calc_dg_cycle.py

**Examples:**

::

  calc_dg_cycle.py -d a_b/out_free -s + b_c/out_free -s + a_c/out_free -s - --dualtopology
  calc_dg_cycle.py -d a_b/out_free -s + b_c/out_free -s + a_c/out_free -s - --singletopology comb

**Description:**

Calculates thermodynamic cycle closure for a set of simulations. This can be performed either for dual topology results, or single topology results. With single topology simulations, the electrostatic and van der Waals results can either be considered separately ``--singletopology sep`` or together, ``--singletopology comb``. 

-----------------------
calc_gcap_surface.py
-----------------------
.. argparse::
   :module: calc_gcap_surface
   :func: get_arg_parser
   :prog: calc_gcap_surface.py

**Examples:**

::

  calc_gcap_surface.py -d out_gcap -v 300. 
  calc_gcap_surface.py -d out_gcap --save-figures -v 300.
  calc_gcap_surface --subdir b_-9.700 --estimators mbar -v 300.

**Description:**

Calculates the free energy from a surface-GCAP simulation. The volume of the GCMC region must be given using the ``-v`` flag. To calculate the free energy at a single B value, use the ``--subdir`` flag with calc_dg.py, and the energy can be calculated with any one dimensional free energy method.

-----------------------
calc_gci.py
-----------------------
.. argparse::
   :module: calc_gci
   :func: get_arg_parser
   :prog: calc_gci.py

**Examples:**

::

  calc_gci.py -d out_gcmc/ -v 130.
  calc_gci.py -d out_gcmc/ -v 130. -l 0.2
  calc_gci.py -d out_gcmc/ -v 130. --save-figures
  calc_gci.py -d out_gcmc/ -v 130. --pin_min


**Description:**

Collection of tools to analyse and visualise GCMC titration data of water using grand canonical integration (GCI). Used to plot average number of waters for a given Adams value, i.e. GCMC titration data, calculate transfer free energies from ideal gas, calculate absolute and relative binding free energies of water, calculate and/or estimate optimal number of bound waters. As described in Ross et al., J. Am. Chem. Soc., 2015, 137 (47), pp 14930-14943. 

Error estimates of free energies and optimal number of waters are based on automatic repeated fitting of the ANN from different random initial parameters. This can be increased with ``--nfits``.

-----------------------
calc_replicapath.py
-----------------------
.. argparse::
   :module: calc_replicapath
   :func: get_arg_parser
   :prog: calc_replicapath.py

**Examples:**

::

  calc_replicapath.py -f out_free/lam-0.*/results -p 0.000 1.000
  calc_replicapath.py -f out_free/lam-0.*/results -p 0.000 0.500 1.000 -o replica_paths.png
  calc_replicapath.py -f out_free/t-*/lam-0.000/results -p 25.0 35.0 45.0 -k temperature

**Description:**

This tools plots the path of different replicas in a replica exchange simulation as a function of simulation time.

If the kind of replicas is from :math:`\lambda` replica exchange the ``replica1`` and ``replica2`` etc should be individual :math:`\lambda`-values to plot. 

If the kind of replicas is from REST or temperature replica exchange the ``replica1`` and ``replica2`` etc should be individual temperatures to plot. 

-----------------------
calc_rmsd.py
-----------------------
.. argparse::
   :module: calc_rmsd
   :func: get_arg_parser
   :prog: calc_rmsd.py

**Examples:**

::

  calc_rmsd.py -i benzene.pdb -f out_bnd/all.pdb -r bnz
  calc_rmsd.py -i benzene.pdb -f out_bnd/all.pdb -r bnz -a c4

**Description:**

This tool calculate the RMSD of a ligand in a simulation.

If the ``atom`` name is given, the tool will calculate the RMSD of that atom with respect to its position in ``pdbfile``. Otherwise, the program will calculate the RMSD of the geometric centre with respect to ``pdbfile``.

A force constant to keep the ligand restrained for free energy calculations is estimated from the RMSD using the equipartition theorem.

-----------------------
calc_series.py
-----------------------
.. argparse::
   :module: calc_series
   :func: get_arg_parser
   :prog: calc_series.py

The tool will estimate the number of independent samples for a given observable in the production part using the method of statistical inefficiency. The equilibration time will also be estimated from a method that maximizes the number uncorrelated samples as suggested on alchemistry.org.

Apart from the raw series, the tool can also plot the running average if the ``--average`` flag is set or the moving average if the ``--moving`` flag is used.

Typically only a single ProtoMS results file will be analysed and plotted. However, for the series ``grad`` and ``agrad`` (the gradient and analytical gradient, respectively), multiple results file can be given. In this case, the gradients for each results file is used to estimate the free energy using thermodynamic integration.


-----------------------
calc_ti_decomposed.py
-----------------------
.. argparse::
   :module: calc_ti_decomposed
   :func: get_arg_parser
   :prog: calc_ti_decomposed.py

**Examples:**

::

  calc_ti.py -d out_free/
  calc_ti.py -d out_free/ -l 0.1 -u 0.9
  calc_ti.py -b out_bnd/ -d out_free --dualtopology
  calc_ti.py -d out_free -g out_gas

**Description:**

This tool calculates free energies of individual energetic components using the method of thermodynamic integration (TI).

The program expects that in the ``directory`` there exist an output folder for each :math:`\lambda`-value, eg. ``lam-0.000`` and ``lam-1.000``


Block estimates can be constructed by combining ``-l`` and ``-u``. For instance, these commands calculates the free energy while incrementally increasing the equilibration ::

  for X in `seq 0.0 0.1 1.0`
  do
  calc_ti_decomposed.py -d out_free -l $x
  done

-----------------------
clear_gcmcbox.py
-----------------------
.. argparse::
   :module: clear_gcmcbox
   :func: get_arg_parser
   :prog: clear_gcmcbox.py

**Examples:**

::

  clear_gcmcbox.py -b gcmc_box.pdb -s water.pdb
  clear_gcmcbox.py -b gcmc_box.pdb -s water.pdb -o water_cleared.pdb

**Description:**

This tool clears a GCMC or JAWS-1 simulation box from any bulk water placed there by the solvation method.

In a GCMC and JAWS-1 simulation the bulk water is prevented to enter or exit a GCMC or JAWS-1 simulation box. Therefore, bulk water that are within this box needs to be removed prior to the GCMC or JAWS-1 simulation. 

The ``boxfile`` is typically created by ``make_gcmcbox.py`` and the ``waterfile`` is typically created by ``solvate.py`` and can be either a droplet or a box.

-----------------------
convertatomnames.py
-----------------------
.. argparse::
   :module: convertatomnames
   :func: get_arg_parser
   :prog: convertatomnames.py

**Examples:**

::

  convertatomnames.py -p protein.pdb
  convertatomnames.py -p protein.pdb -c $PROTOMSHOME/data/atomnamesmap.dat
  convertatomnames.py -p protein.pdb -s charmm

**Description:**

This tool converts residue and atom names to ProtoMS convention. 

This script modfies in particular names of hydrogen atoms, but also some residue names, e.g. histidines.

A file containing conversion instructions for amber and charmm is available in the ``$PROTOMSHOME/data`` folder.


-----------------------
convertwater.py
-----------------------
.. argparse::
   :module: convertwater
   :func: get_arg_parser
   :prog: convertwater.py

**Examples:**

::

  convertwater.py -p protein.pdb
  convertwater.py -p protein.pdb -m tip3p
  convertwater.py -p protein.pdb --ignoreh

**Description:**

This tool converts water molecules to a specific model.

Currently the script recognizes TIP3P and TIP4P water models. The valid values for ``style`` is therefore ``t4p, tip4p, tp4, t3p, tip3p, tp3``

If the ``--ignoreh`` flag is given, the script will discard the hydrogen atoms found in ``pdbfile`` and add them at a random orientation.


-----------------------
distribute_waters.py
-----------------------
.. argparse::
   :module: distribute_waters
   :func: get_arg_parser
   :prog: distribute_waters.py

**Examples:**

::

  distribute_waters.py -b 53.4 56.28 13.23 10 10 10 -m 12
  distribute_waters.py -b 53.4 56.28 13.23 10 10 10 -m 12 --model t3p --resname T3P
  distribute_waters.py -b 53.4 56.28 13.23 10 10 10 -m myonewater.pdb --number 12 -o mywatersinbox.pdb


**Description:**

This tool can place water molecules at random within a GCMC or JAWS-1 simulation box.

It can place molecules in random positions and orientations with their geometry center restricted to the given dimensions of a box.


-----------------------
divide_pdb.py
-----------------------
.. argparse::
   :module: divide_pdb
   :func: get_arg_parser
   :prog: divide_pdb.py

**Examples:**

::
  divide_pdb.py
  divide_pdb.py -i mypmsout.pdb -o individual -p outfolder/ 


**Description:**

This tool splits up a PDB file with multiple models (the keyword END defines the end of a model) into several PDB files.


-----------------------
generate_input.py
-----------------------
.. argparse::
   :module: generate_input
   :func: get_arg_parser
   :prog: generate_input.py

**Examples:**

::

  generate_input.py -s dualtopology -l lig1.pdb lig2.pdb -p protein.pdb -t li1-li2.tem -pw droplet.pdb -lw lig1_wat.pdb --lambas 8
  generate_input.py -s dualtopology -l lig1.pdb dummy.pdb -t li1-dummy.tem -lw lig1_wat.pdb --absolute
  generate_input.py -s gcmc -p protein.pdb -pw droplet.pdb --adams -4 -2 0 2 4 6 --gcmcwater gcmc_water.pdb --gcmcbox gcmc_box.pdb
  generate_input.py -s sampling -l lig1.pdb -t lig1.tem --dovacuum

**Description:**

This tool generates input files with commands for ProtoMS.

The settings generate are made according to experience and should work in most situations.

The tool will create at most two ProtoMS command files, one for the protein simulation and one for the ligand simulation. These can be used to run ProtoMS, e.g. ::

  $PROTOMS/protoms3 run_free.cmd

-----------------------
make_dummy.py
-----------------------
.. argparse::
   :module: make_dummy
   :func: get_arg_parser
   :prog: make_dummy.py

**Examples:**

::

  make_dummy.py -f benzene.pdb
  make_dummy.py -f benzene.pdb -o benzene_dummy.pdb

**Description:**

This tool makes a matching dummy particle for a solute.

The dummy particle will be placed at the centre of the solute.


-----------------------
make_gcmcbox.py
-----------------------
.. argparse::
   :module: make_gcmcbox
   :func: get_arg_parser
   :prog: make_gcmcbox.py

**Examples:**

::

  make_gcmcbox.py -s benzene.pdb
  make_gcmcbox.py -s benzene.pdb -p 0.0
  make_gcmcbox.py -s benzene.pdb -o benzene_gcmc_box.pdb

**Description:**

This tool makes a GCMC or JAWS-1 simulation box to fit on top of a solute.

The box will be created so that it has the extreme dimensions of the solute and then ``padding`` will be added in each dimension

The box can be visualised with most common programs, e.g. ::

  vmd -m benzene.pdb benzene_gcmc_box.pdb

this is a good way to see that the box is of appropriate dimensions.

When an appropriate box has been made, it can be used by ``solvate.py`` to fill it with water.

-----------------------
make_single.py
-----------------------
.. argparse::
   :module: make_single
   :func: get_arg_parser
   :prog: make_single.py

**Examples:**

::

  make_single.py -t0 benzene.tem -t1 toluene.tem -p0 benzene.pdb -p1 toluene.pdb
  make_single.py -t0 benzene.tem -t1 toluene.tem -p0 benzene.pdb -p1 toluene.pdb -m bnz2tol.dat
  make_single.py -t0 benzene.tem -t1 toluene.tem -p0 benzene.pdb -p1 toluene.pdb -o bnz-tol

**Description:**


This tool makes ProtoMS template files for single topology free energy simulations.

The program will automatically try to match atoms in ``template0`` with atoms in ``template1``. It will do this by looking for atoms with the same atom type that are on top of each other in ``pdbfile0`` and ``pdbfile1``. A cut-off of 0.02 A2 will be used for this. All atoms that cannot be identified in this way are written to the screen and the user has to enter the corresponding atoms. If no corresponding atom exists, i.e., the atom should be perturbed to a dummy, the user may enter blank. 

The user may also write the corresponding atoms to a file and provide it as ``map`` above. In this file there should be one atom pair on each line, separated by white-space. A dummy atom should be denoted as ``DUM``. If ``map`` is not given, the program will write the created correspondence map to a file based on the ``outfile`` string.

Currently, dummy atoms are not supported in the solute at :math:`\lambda=0.0`. Therefore, this solute needs to be the larger one.

The tool will write two ProtoMS template files, one for the electrostatic perturbation, one for the van der Waals perturbation and one for the combined perturbation. These template files will end in ``_ele.tem``, ``_vdw.tem``, ``_comb.tem`` respectively. 

A summary of the charges and van der Waals parameters in the four states will be printed to the screen. This information should be checked carefully. 


-----------------------
merge_templates.py
-----------------------
.. argparse::
   :module: merge_templates
   :func: get_arg_parser
   :prog: make_templates.py

**Examples:**

::

  merge_templates.py -f benzene.tem dummy.tem -o bnz-dummy.tem


**Description:**

This tool combines several ProtoMS template files into a single template file.

The force field parameters in ``file2`` will be re-numbered so that they do not conflict with ``file1``. This is important when you want to load both parameters into ProtoMS at the same time.

-----------------------
plot_theta.py
-----------------------
.. argparse::
   :module: plot_theta
   :func: get_arg_parser
   :prog: plot_theta.py
	  

**Examples:**

::

  plot_theta.py -m WA1 --skip 50
  plot_theta.py -m WA1 -p theta_wa1


**Description:**

This tool plots the theta distribution resulting from a JAWS stage one simulation.

Two different histograms will be generated. One in which all different copies of the same molecule are added up, and a different one where each copy is displayed individually.

-----------------------
scoop.py
-----------------------
.. argparse::
   :module: scoop
   :func: get_arg_parser
   :prog: scoop.py

**Examples:**

::

  scoop.py -p protein.pdb
  scoop.py -p protein.pdb  -l benzene.pdb
  scoop.py -p protein.pdb  --center "0.0 0.0 0.0"
  scoop.py -p protein.pdb  --center origin.dat
  scoop.py -p protein.pdb  --innercut 10 --outercut 16
  scoop.py -p protein.pdb  --exclude 189 190
  scoop.py -p protein.pdb  --added 57 58 59 


**Description:**

This tool truncates a protein and thereby creating a scoop.

All residues outside ``ocut`` is removed completely. ``icut`` is used to separate the scoop model into two different regions, that possibly can have different sampling regimes. The sampling regimes are determined by ``--flexin`` and ``--flexout``. 

If the user would like to finetune the residues in the scoop this can be done with ``--excluded`` to discard specific residues or ``--added`` to include specific residues.

The scoop will be centred on the ``ligandfile`` is such a file is provided. Otherwise, it will be centred on the flag ``--center``. The argument to this flag can be either a string with three numbers specifying the centre, as in example three above. It can also be the name of a file containing the centre, as in example four above.

Crystallographic waters that are in ``proteinfile`` will also be truncated at ``ocut``

The PDB file will contain specific instructions for ProtoMS to automatically enforce the values of  ``--flexin`` and ``--flexout``.



-----------------------
solvate.py
-----------------------
.. argparse::
   :module: solvate
   :func: get_arg_parser
   :prog: solvate.py

**Examples:**

::

  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -s benzene.pdb
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -s benzene.pdb -p 12.0 
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -s benzene.pdb -pr protein.pdb -g droplet
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -s benzene.pdb -pr protein.pdb -g droplet -r 24.0
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -pr protein.pdb -g droplet -c 0.0
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -pr protein.pdb -g droplet -c "0.0 10.0 20.0"
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -pr protein.pdb -g droplet -c "76 86"
  solvate.py -b $PROTOMSHOME/data/wbox_tip4p.pdb -s gcmc_box.pdb -g flood


**Description:**

This tool solvates a ligand in either a droplet or a box of water. It can also flood a GCMC or JAWS-1 simulations box with waters.

Pre-equilibrated boxes to use can be found in the ``$PROTOMSHOME/data`` folder.

To solvate small molecule it is sufficient to give the ``solutefile`` as in the first example above. This produces a box with at least 10 A between the solute and the edge of the water box, which should be sufficient in most situation. Use ``padding`` to increase or decrease the box size as in the second example. The solvation box is created by replicating the pre-equilibrated box in all dimensions and then removing waters that overlap with solute atoms.

To solvate a protein in a droplet, specify ``proteinfile`` and ``droplet`` as in the third example above. This produces a droplet with radius of 30 A, which was chosen to work well with the default options in ``scoop.py``. Use ``radius`` to obtain a smaller or larger droplet as in the fourth example. The centre of the droplet can be on a ligand if ``ligandfile`` is specified. Otherwise, the ``center``argument is used. This argument can be either ``cent`` (the default) that places the droplet at the centre of the protein. It can also take a single number as in the fifth example above in case it is placed at this coordinate in all dimensions. It can also take a string with three numbers which is the origin of the droplet in x, y, and z dimensions, see the sixth example above. If two numbers are given as in the seventh example above, it is assumed that this is an atom range and the droplet will be placed at the centre of these atoms. The droplet is created by putting random waters from the pre-equilibrated box on a grid, displacing them slightly in a random fashion.

The tool can also be used to fill a box with waters for GCMC and JAWS-1 simulations, similar to ``distribute_waters.py``. In this case the solute is typically a box created by ``make_gcmcbox.py`` and ``flood`` needs to be specified, see the last example above. This gives a box filled with the bulk number of waters.


-----------------------
split_jawswater.py
-----------------------
.. argparse::
   :module: split_jawswater
   :func: get_arg_parser
   :prog: split_jawswater.py

**Examples:**

::

  split_jawswater.py -w waters.pdb
  split_jawswater.py -w waters.pdb -o jaws2_


**Description:**

This tool splits a PDB file containing multiple water molecules into PDB files appropriate for JAWS-2. 

For each water molecule in ``pdbfile`` the tool will write a PDB file with individual water molecules named ``outprefix+watN.pdb`` where N is the serial number of the water molecule. Furthermore, the tool will write a PDB file with all the other molecules and name if ``outprefix+notN.pdb`` where again N is the serial number of the water molecule. In these latter PDB-files, the water residue name is changed to that of the bulk water, e.g., ``t3p`` or ``t4p``.

For instance, if ``waters.pdb`` in the second example above contains 3 water molecule, this tool will create the following files: ::

  jaws2_wat1.pdb
  jaws2_wat2.pdb
  jaws2_wat3.pdb

  jaws2_not1.pdb
  jaws2_not2.pdb
  jaws2_not3.pdb

