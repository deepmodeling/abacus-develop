if(ENABLE_LCAO)
    list(APPEND objects
        evolve_elec.cpp
        evolve_psi.cpp
        band_energy.cpp
        middle_hamilt.cpp
        norm_psi.cpp
        propagator.cpp
        propagator_cn2.cpp
        propagator_taylor.cpp
        propagator_etrs.cpp
        upsi.cpp
        td_velocity.cpp
        td_current.cpp
        snap_psibeta_half_tddft.cpp
    )

    add_library(
        tddft
        OBJECT
        ${objects}
    )

    if(ENABLE_COVERAGE)
    add_coverage(tddft)
    endif()

    IF (BUILD_TESTING)
        if(ENABLE_MPI)
            add_subdirectory(test)
        endif()
    endif()

endif()
