add_subdirectory(module_rt)
add_subdirectory(module_deepks)
add_subdirectory(module_dftu)
add_subdirectory(module_hcontainer)
add_subdirectory(module_deltaspin)

if(ENABLE_LCAO)
    add_subdirectory(module_operator_lcao)
    list(APPEND objects
        hamilt_lcao.cpp
        module_operator_lcao/operator_lcao.cpp
        module_operator_lcao/veff_lcao.cpp
        module_operator_lcao/meta_lcao.cpp
        module_operator_lcao/op_dftu_lcao.cpp
        module_operator_lcao/deepks_lcao.cpp
        module_operator_lcao/op_exx_lcao.cpp
        module_operator_lcao/overlap_new.cpp
        module_operator_lcao/ekinetic_new.cpp
        module_operator_lcao/nonlocal_new.cpp
        module_operator_lcao/td_ekinetic_lcao.cpp
        module_operator_lcao/td_nonlocal_lcao.cpp
        module_operator_lcao/td_pot_hybrid.cpp
        module_operator_lcao/dspin_lcao.cpp
        module_operator_lcao/dftu_lcao.cpp
        pulay_force_stress_center2.cpp
        FORCE_STRESS.cpp
        FORCE_gamma.cpp
        FORCE_k.cpp
        stress_tools.cpp
        edm.cpp
        grid_init.cpp
        spar_dh.cpp
        spar_exx.cpp
        spar_hsr.cpp
        spar_st.cpp
        spar_u.cpp
		LCAO_set_fs.cpp
		LCAO_set_st.cpp
		LCAO_nl_mu.cpp	
        LCAO_set_zero.cpp
        LCAO_allocate.cpp
        LCAO_set_mat2d.cpp
        LCAO_init_basis.cpp
        setup_exx.cpp
        setup_deepks.cpp
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
        ${objects}
    )

    if(ENABLE_COVERAGE)
    add_coverage(hamilt_lcao)
    endif()

endif()


