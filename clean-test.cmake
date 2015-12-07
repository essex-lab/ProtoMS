set(DELFILES
    # setup test
    "test_setup/dcb.prepi"
    "test_setup/dcb.frcmod"
    "test_setup/dcb.zmat"
    "test_setup/dcb.tem"
    "test_setup/dcb_box.pdb"
    "test_setup/protein_scoop.pdb"
    "test_setup/water.pdb"
    # equilibration test
    "test_equil/run_bnd.cmd"
    "test_equil/out_bnd"
    # sampling test
    "test_sampling/run_bnd.cmd"
    "test_sampling/out_bnd"
    )

foreach(f ${DELFILES})
    set(f ../tests/${f})
    if(EXISTS ${f})
        message("Removing ${f}")
        file(REMOVE_RECURSE ${f})
    endif()
endforeach()
