*************
Test Suite
*************

The ProtoMS test suite can be found in the ``$PROTOMSHOME/tests`` directory. It contains a set of Python scripts and all required input files and reference output files to run a sanity check on the ProtoMS code, both the source (Fortran) code and the (Python) tools. In this page you will find a list of the different tests, a brief indication of which part of ProtoMS each of the tests is checking and instructions to run each of the individual tests separately, or all of them as a whole.

==========================================
Dependencies
==========================================

The Python module ``nose`` is required to run the test suite. You can find more information on nose on its website ``nose.readthedocs.org/en/latest/``.

==========================================
Running all tests
==========================================

The simplest and recommended way to run the tests is to run the command ``ctest`` while in the build directory ``$PROTOMSHOME/build``.  This will run all tests and report the success or failure of each.  For more information use the command ``ctest -V`` which will print all output from the tests as well as output from both the Python and Fortran components of ProtoMS.

If ProtoMS was compiled without MPI, the following test scripts will not be run:

* ``test_mpi_install.py``
* ``test_gcmc.py``
* ``test_jaws2_sim.py``
* ``test_reti_sngl.py``
* ``test_reti_dbl.py``

==========================================
Individual tests
==========================================

In this section you will a list of all tests, with a brief explanaition of which part of the ProtoMS code they are testing.

Most of the test scripts define multiple individual tests, a setup test and a simulation test.
From the ProtoMS build directory (i.e. ``$PROTOMSHOME/build``), to run all tests within a single test script use ``python ../tests/test_<scriptname>.py``.
To run only a single test stage use ``python ../tests/test_<scriptname>.py <FullTestName>`` or ``ctest -R <FullTestName>``

----------------------------
List of tests
----------------------------

A list of all Python scripts that correspond to each of the tests is shown below:

* ``test_install_dependencies.py``
* ``test_ligand_setup.py``
* ``test_parameters_ff.py``
* ``test_path.py``
* ``test_tools_protoms.py``
* ``test_prot_setup.py``
* ``test_equil.py``
* ``test_energies.py``
* ``test_sampling.py``
* ``test_jaws1.py``
* ``test_mpi_install.py``
* ``test_gcmc.py``
* ``test_jaws2.py``
* ``test_reti_sngl.py``
* ``test_reti_dbl.py``


----------------------------
test_install_dependencies.py
----------------------------

**Coverage**

  This test covers the requirements for the installation of ProtoMS.  It checks that the AmberTools are installed and available at ``$AMBERHOME`` and that the Python modules numpy, scipy and matplotlib are available.

----------------------------
test_ligand_setup.py
----------------------------

**Coverage**

  Checks that the set up of ligands with the ProtoMS tools generates the expected results.

**Reference Data Location**

 ``$PROTOMSHOME/tests/setup``

----------------------------
test_parameters_ff.py
----------------------------

**Coverage**

  Checks that all expected parameter files are found in ``$PROTOMSHOME/parameter``.

----------------------------
test_path.py
----------------------------

**Coverage**

  Checks that ``$PROTOMSHOME`` has been set correctly.

----------------------------
test_tools_protoms.py
----------------------------

**Coverage**

  Checks that all expected Python scripts corresponding to the ProtoMS tools are present in ``$PROTOMSHOME/tools``.

----------------------------
test_prot_setup.py
----------------------------

**Contains**
* ProtSetupTest

**Coverage**

  Checks that the set up of protein and ligand with the ProtoMS tools generates the expected results.

**Reference Data Location**

 ``$PROTOMSHOME/tests/setup/``

----------------------------
test_equil_prot.py
----------------------------

**Contains**
* EquilSetupTest
* EquilSimulationTest

**Coverage**

  Checks both setup and run of the ``equilibration`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/equil/``

----------------------------
test_energies.py
----------------------------

**Contains**
* EnergiesSimulationTip3pTest
* EnergiesSimulationTip4pTest

**Coverage**

  Checks the generation of correct energies for different water models used as solvent.

**Reference Data Location**

 ``$PROTOMSHOME/tests/energies/``

----------------------------
test_sampling.py
----------------------------

**Contains**
* SamplingSetupTest
* SamplingSimulationTest

**Coverage**

  Checks both setup and run of the ``sampling`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/sampling/``

----------------------------
test_jaws1.py
----------------------------

**Contains**
* Jaws1SetupTest
* Jaws1SimulationTest

**Coverage**

  Checks both setup and run of the ``jaws1`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/jaws1/``

----------------------------
test_mpi_install.py
----------------------------

**Coverage**

  Checks that MPI is available for running simulations requiring it.

----------------------------
test_gcmc.py
----------------------------

**Contains**
GcmcSetupBoxTest
GcmcSetupTest
GcmcSimulationTest

**Coverage**

  Checks both setup and run of the ``gcmc`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/gcmc/``

----------------------------
test_jaws2.py
----------------------------

**Contains**
* Jaws2SetupTest
* Jaws2SimulationTest

**Coverage**

  Checks both setup and run of the ``jaws2`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/jaws2/``

----------------------------
test_reti_sngl.py
----------------------------

**Contains**
* RetiSnglSetupTest
* RetiSnglSimulationFreeTest
* RetiSnglSimulationGasTest

**Coverage**

  Checks both setup and run of the ``singletopology`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/RETI_sngl/``

----------------------------
test_reti_dbl.py
----------------------------

**Contains**
* RetiDblSetupTest
* RetiDblSimulationTest

**Coverage**

  Checks both setup and run of the ``dualtopology`` simulation type among those offered by ``protoms.py``.

**Reference Data Location**

 ``$PROTOMSHOME/tests/RETI_dbl/``
