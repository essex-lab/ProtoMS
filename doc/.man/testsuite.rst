*************
Test Suite
*************

The ProtoMS test suite can be found in the ``$PROTOMSHOME/test``. It contains a set of python scripts and all required input files to run a sanity check on the ProtoMS code, both the source (Fortran) code and the (python) tools. In this page you will find a list of the different available tests, a brief indication of which part of ProtoMS each of the tests is checking, instructions to run each of the individual tests separately, or all of them as a whole.

==========================================
Dependencies
==========================================

Nose is required to run the test suite. You can find more information on nose on its website ``nose.readthedocs.org/en/latest/``.

==========================================
Running all tests
==========================================

All tests can be run at once using CTest. To proceed, after compilation, go to the build directory ``$PROTOMSHOME/buil`` and type the command:
* ``make test``
or, alternatively
* ``ctest``

==========================================
Individual tests
==========================================

In this section you will a list of all tests, with a brief explanaition of which part of the ProtoMS code they are testing and the instructions to run them.

----------------------------
List of tests
----------------------------

A list of all python scripts that correspond to each of the tests is shown below:

* ``test_install_dependencies.py``
* ``test_ligand_setup.py``
* ``test_mpi_install.py``
* ``test_parameters_ff.py``
* ``test_path.py``
* ``test_tools_protoms.py``
* ``test_prot_setup.py``
* ``test_equil_prot.py``
* ``test_energies.py``
* ``test_sampling_prot.py``
* ``test_gcmc_sim.py``
* ``test_jaws1_sim.py``
* ``test_jaws2_sim.py``
* ``test_reti_sngl.py``
* ``test_reti_dbl.py``


----------------------------
test_install_dependencies.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_install_dependencies.py``

**Coverage**

  This test covers the requirements for the installation of ProtoMS, both in terms of library depen

----------------------------
test_ligand_setup.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_ligand_setup.py``

**Coverage**

  This test covers that the set up of ligands with the ProtoMS tools generates the expected results.

**Input**

* ``$PROTOMSHOME/test_setup/dcb.pdb``


----------------------------
test_mpi_install.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_mpi_install.py``

**Coverage**

  This test covers that the requirements for ProtoMS to run with opem_mpi.

----------------------------
test_parameters_ff.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_parameters_ff.py``

**Coverage**

  This test checks that all expected parameter files are found in ``$PROTOMSHOME/parameter``.

----------------------------
test_path.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_path.py``

**Coverage**

  This test checks that ``$PROTOMSHOME`` has been set correctly.

----------------------------
test_tools_protoms.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/``

**Use**

 ``nosetests tests/test_tools_protoms.py``

**Coverage**

  This test checks that all expected python scripts corresponding to the ProtoMS tools are present in ``$PROTOMSHOME/tools``.

----------------------------
test_prot_setup.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_setup/``

**Use**

 ``nosetests tests/test_setup/test_prot_setup.py``

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

**Use**

 ``nosetests tests/test_equil/test_equil_prot.py``

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

**Use**

 ``nosetests tests/test_energies/test_energies.py``

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

**Use**

 ``nosetests tests/test_sampling/test_sampling_prot.py``

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

**Use**

 ``nosetests tests/test_gcmc/test_gcmc_sim.py``

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

**Use**

 ``nosetests tests/test_jaws1/test_jaws1_sim.py``

**Coverage**

  This test checks both setup and run of the ``jaws1`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_jaws1/protein.pdb``
* ``$PROTOMSHOME/test_jaws1/fragment.pdb``
* ``$PROTOMSHOME/test_jaws1/water.pdb``

----------------------------
test_jaws2_sim.py
----------------------------

**Location**

 ``$PROTOMSHOME/tests/test_jaws2/``

**Use**

 ``nosetests tests/test_jaws2/test_jaws2_sim.py``

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

**Use**

 ``nosetests tests/test_RETI_sngl/test_reti_sngl.py``

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

**Use**

 ``nosetests tests/test_RETI_dbl/test_reti_dbl.py``

**Coverage**

  This test checks both setup and run of the ``dualtopology`` simulation type among those offered by ``protoms.py``.

**Input**

* ``$PROTOMSHOME/test_RETI_dbl/ethane.pdb``
* ``$PROTOMSHOME/test_RETI_dbl/methanol.pdb``




