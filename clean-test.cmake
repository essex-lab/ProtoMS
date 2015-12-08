set(DELFILES
    # Setup test
    "test_setup/dcb.prepi"
    "test_setup/dcb.frcmod"
    "test_setup/dcb.zmat"
    "test_setup/dcb.tem"
    "test_setup/dcb_box.pdb"
    "test_setup/protein_scoop.pdb"
    "test_setup/water.pdb"
    # Equilibration test
    "test_equil/run_bnd.cmd"
    "test_equil/out_bnd"
    # Sampling test
    "test_sampling/run_bnd.cmd"
    "test_sampling/out_bnd"
    # GCMC test
    "test_gcmc/gcmc_box.pdb"
    "test_gcmc/gcmc_wat.pdb"
    "test_gcmc/run_bnd.cmd"
    "test_gcmc/out_gcmc"
    # JAWS1 test
    "test_jaws1/water_clr.pdb"
    "test_jaws1/jaws1_box.pdb"
    "test_jaws1/jaws1_wat.pdb"
    "test_jaws1/run_jaws.cmd"
    "test_jaws1/out_jaws"
    )

foreach(f ${DELFILES})
    set(f ../tests/${f})
    if(EXISTS ${f})
        message("Removing ${f}")
        file(REMOVE_RECURSE ${f})
    endif()
endforeach()
