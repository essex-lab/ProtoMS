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
    # JAWS2 test
    "test_jaws2/jaws2_wat3.pdb"
    "test_jaws2/jaws2_wat2.pdb"
    "test_jaws2/jaws2_wat1.pdb"
    "test_jaws2/jaws2_not3.pdb"
    "test_jaws2/jaws2_not2.pdb"
    "test_jaws2/jaws2_not1.pdb"
    "test_jaws2/jaws1_box.pdb"
    "test_jaws2/jaws1_wat.pdb"
    "test_jaws2/run_jaws2-w1_bnd.cmd"
    "test_jaws2/run_jaws2-w2_bnd.cmd"
    "test_jaws2/run_jaws2-w3_bnd.cmd"
    "test_jaws2/out_jaws2-w1"
    # RETI single test
    "test_RETI_sngl/eth-meo_comp.tem"
    "test_RETI_sngl/eth-meo_ele.tem"
    "test_RETI_sngl/eth-meo_gas.tem"
    "test_RETI_sngl/run_free.cmd"
    "test_RETI_sngl/run_gas.cmd"
    "test_RETI_sngl/out_free"
    "test_RETI_sngl/out_gas"
    # RETI double test
    "test_RETI_dbl/eth-meo.tem"
    "test_RETI_dbl/run_free.cmd"
    "test_RETI_dbl/out_free"
    )

foreach(f ${DELFILES})
    set(f ../tests/${f})
    if(EXISTS ${f})
        message("Removing ${f}")
        file(REMOVE_RECURSE ${f})
    endif()
endforeach()
