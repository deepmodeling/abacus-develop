if (ENABLE_LIBRI)

    list(APPEND objects
        irreducible_sector.cpp
        irreducible_sector_bvk.cpp
        symmetry_rotation.cpp
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