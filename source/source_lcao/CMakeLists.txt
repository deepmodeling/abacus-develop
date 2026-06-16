add_subdirectory(module_rt)
add_subdirectory(module_deepks)
add_subdirectory(module_dftu)
add_subdirectory(module_deltaspin)

if(ENABLE_LCAO)
    add_subdirectory(module_operator_lcao)
    add_subdirectory(module_extrap)
    set(hamilt_lcao_sources
        hamilt_lcao.cpp
        dftu_lcao.cpp
        pulay_fs_center2.cpp
        FORCE_STRESS.cpp
        FORCE_gamma.cpp
        FORCE_k.cpp
        stress_tools.cpp
        edm.cpp
        spar_dh.cpp
        spar_exx.cpp
        spar_hsr.cpp
        spar_st.cpp
        spar_u.cpp
        LCAO_set.cpp
        LCAO_set_fs.cpp
        LCAO_set_st.cpp
        LCAO_nl_mu.cpp
        LCAO_set_zero.cpp
        LCAO_allocate.cpp
        LCAO_set_mat2d.cpp
        LCAO_init_basis.cpp
        setup_nonlocal.cpp
        setup_exx.cpp
        setup_deepks.cpp
        setup_dm.cpp
        rho_tau_lcao.cpp
        record_adj.cpp
        center2_orb.cpp
        center2_orb-orb11.cpp
        center2_orb-orb21.cpp
        center2_orb-orb22.cpp
        wavefunc_in_pw.cpp
    )

    add_library(
        hamilt_lcao
        OBJECT
        ${hamilt_lcao_sources}
    )

    if(ENABLE_COVERAGE)
    add_coverage(hamilt_lcao)
    endif()

    IF (BUILD_TESTING)
        add_subdirectory(test)
    endif()

endif()


