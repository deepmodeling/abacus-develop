add_subdirectory(operator_lcao)

if(ENABLE_LCAO)
    list(APPEND objects
        hamilt_lcao.cpp
        operator_lcao/operator_lcao.cpp
        operator_lcao/ekinetic_lcao.cpp
        operator_lcao/nonlocal_lcao.cpp
        operator_lcao/overlap_lcao.cpp
        operator_lcao/veff_lcao.cpp
        operator_lcao/meta_lcao.cpp
        operator_lcao/op_dftu_lcao.cpp
        operator_lcao/deepks_lcao.cpp
        operator_lcao/op_exx_lcao.cpp
        operator_lcao/lambda_lcao.cpp
        FORCE_STRESS.cpp
        FORCE_gamma.cpp
        FORCE_gamma_edm.cpp
        FORCE_gamma_tvnl.cpp
        FORCE_gamma_vl.cpp
        FORCE_k.cpp
        LCAO_gen_fixedH.cpp
        LCAO_hamilt.cpp
        LCAO_matrix.cpp
        LCAO_nnr.cpp
        record_adj.cpp
        center2_orb-orb11.cpp
        center2_orb-orb21.cpp
        center2_orb-orb22.cpp
        DM_gamma.cpp
        DM_gamma_2d_to_grid.cpp
        DM_k.cpp
        local_orbital_charge.cpp
        local_orbital_wfc.cpp
        dm_2d.cpp
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