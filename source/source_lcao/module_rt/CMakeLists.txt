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
        td_info.cpp
        velocity_op.cpp
        snap_psibeta_half_tddft.cpp
        td_folding.cpp
        solve_propagation.cpp
        boundary_fix.cpp
    )

    if(USE_CUDA)
        list(APPEND objects
            kernels/cuda/snap_psibeta_kernel.cu
            kernels/cuda/snap_psibeta_gpu.cu
            kernels/cuda/norm_psi_kernel.cu
            kernels/cuda/band_energy_kernel.cu
        )
    endif()

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
