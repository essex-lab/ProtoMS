*************
Test Suite
*************

The ProtoMS test suite can be found in the ``$PROTOMSHOME/tests`` directory. It contains a set of python scripts and all required input files to run a sanity check on the ProtoMS code, both the source (Fortran) code and the (python) tools. In this page you will find a list of the different tests, a brief indication of which part of ProtoMS each of the tests is checking and instructions to run each of the individual tests separately, or all of them as a whole.

==========================================
Dependencies
==========================================

The python module ``nose`` is required to run the test suite. You can find more information on nose on its website ``nose.readthedocs.org/en/latest/``.

==========================================
Running all tests
==========================================

The simplest and recommended way to run the tests is to run the command ``ctest`` while in the build directory ``$PROTOMSHOME/build``.  This will run all tests and report the success or failure of each.  For more information use the command ``ctest -V`` which will print all output from the tests as well as output from both the python and Fortran components of ProtoMS.  To clean up after running the tests use the command ``make clean-test`` from the build directory.

If ProtoMS was compiled without MPI, the following tests will not be run automatically:

* ``test_mpi_install.py``
* ``test_jaws2_sim.py``
* ``test_reti_sngl.py``
* ``test_reti_dbl.py``

==========================================
Individual tests
==========================================

In this section you will a list of all tests, with a brief explanaition of which part of the ProtoMS code they are testing.  Individual tests may be run either by using the command ``ctest -R test_testname`` from the build directory, or by running the test script directly using ``python test_testname.py`` in the directories indicated below e.g. ``$PROTOMSHOME/tests/test_gcmc``.

----------------------------
List of tests
----------------------------

A list of all python scripts that correspond to each of the tests is shown below:

* ``test_install_dependencies.py``
* ``test_ligand_setup.py``
* ``test_parameters_ff.py``
* ``test_path.py``
* ``test_tools_protoms.py``
* ``test_prot_setup.py``
* ``test_equil_prot.py``
* ``test_energies.py``
* ``test_sampling_prot.py``
* ``test_gcmc_sim.py``
* ``test_jaws1_sim.py``
* ``test_mpi_install.py``
* ``test_jaws2_sim.py``
* ``test_reti_sngl.py``
* ``test_reti_dbl.py``


----------------------------
test_install_dependencies.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test covers the requirements for the installation of ProtoMS.  It checks that the AmberTools are installed and available at ``$AMBERHOME`` and that the python modules numpy, scipy and matplotlib are available.

----------------------------
test_ligand_setup.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test checks that the set up of ligands with the ProtoMS tools generates the expected results.

**Input**

* ``$PROTOMSHOME/test_setup/dcb.pdb``

----------------------------
test_parameters_ff.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test checks that all expected parameter files are found in ``$PROTOMSHOME/parameter``.

----------------------------
test_path.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test checks that ``$PROTOMSHOME`` has been set correctly.

----------------------------
test_tools_protoms.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test checks that all expected python scripts corresponding to the ProtoMS tools are present in ``$PROTOMSHOME/tools``.

----------------------------
test_prot_setup.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_setup/``

**Coverage**

  This test checks that the set up of protein and ligand with the ProtoMS tools generates the expected results.

**Input**

* ``$PROTOMSHOME/test_setup/dcb.pdb``
* ``$PROTOMSHOME/test_setup/protein.pdb``

----------------------------
test_equil_prot.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_equil/``

**Coverage**

  This test checks both setup and run of the ``equilibration`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_equil/dcb.pdb``
* ``$PROTOMSHOME/test_equil/protein.pdb``

----------------------------
test_energies.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_energies/``

**Coverage**

  This test checks checks the generation of the correct energies for different water models used as solvent.

**Input**

* ``$PROTOMSHOME/test_test_energies/t3p.pdb``
* ``$PROTOMSHOME/test_test_energies/t4p.pdb``
* ``$PROTOMSHOME/test_test_energies/run_t3p.cmd``
* ``$PROTOMSHOME/test_test_energies/run_t4p.cmd``
* ``$PROTOMSHOME/test_setup/protein_scoop.pdb``

----------------------------
test_sampling_prot.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_sampling/``

**Coverage**

  This test checks both setup and run of the ``sampling`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_sampling/dcb.pdb``
* ``$PROTOMSHOME/test_sampling/protein.pdb``

----------------------------
test_gcmc_sim.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_gcmc/``

**Coverage**

  This test checks both setup and run of the ``gcmc`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_gcmc/protein.pdb``
* ``$PROTOMSHOME/test_gcmc/wat.pdb``
* ``$PROTOMSHOME/test_gcmc/gcmc_box.pdb``
* ``$PROTOMSHOME/test_gcmc/water.pdb``

----------------------------
test_jaws1_sim.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_jaws1/``

**Coverage**

  This test checks both setup and run of the ``jaws1`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_jaws1/protein.pdb``
* ``$PROTOMSHOME/test_jaws1/fragment.pdb``
* ``$PROTOMSHOME/test_jaws1/water.pdb``

----------------------------
test_mpi_install.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Coverage**

  This test checks that MPI is available for running simulations requiring it.

----------------------------
test_jaws2_sim.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_jaws2/``

**Coverage**

  This test checks both setup and run of the ``jaws2`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_jaws2/protein.pdb``
* ``$PROTOMSHOME/test_jaws2/fragment.pdb``
* ``$PROTOMSHOME/test_jaws2/water.pdb``
* ``$PROTOMSHOME/test_jaws2/jaws2_waters.pdb``

----------------------------
test_reti_sngl.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_RETI_sngl/``

**Coverage**

  This test checks both setup and run of the ``singletopology`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_RETI_sngl/ethane.pdb``
* ``$PROTOMSHOME/test_RETI_sngl/methanol.pdb``
* ``$PROTOMSHOME/test_RETI_sngl/single_cmap.dat``

----------------------------
test_reti_dbl.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_RETI_dbl/``

**Coverage**

  This test checks both setup and run of the ``dualtopology`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_RETI_dbl/ethane.pdb``
* ``$PROTOMSHOME/test_RETI_dbl/methanol.pdb``




