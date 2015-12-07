set(DELFILES
    "test_equil/run_bnd.cmd"
    "test_equil/out_bnd")

foreach(f ${DELFILES})
    set(f ../tests/${f})
    if(EXISTS ${f})
        message("Removing ${f}")
        file(REMOVE_RECURSE ${f})
    endif()
endforeach()
