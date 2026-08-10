if (ENABLE_LIBRI)

    list(APPEND objects
        irreducible_sector.cpp
        irred_sec_bvk.cpp
        symmetry_rotation.cpp
        symm_rot_out.cpp
    )
    add_library(
        module_exx_symmetry
        OBJECT
        ${objects}
    )
    if(BUILD_TESTING)
        add_subdirectory(test)
    endif()

    if(ENABLE_COVERAGE)
        add_coverage(module_exx_symmetry)
    endif()

  endif()