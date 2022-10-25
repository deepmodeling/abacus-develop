add_subdirectory(ks_pw)
add_subdirectory(ks_lcao)
list(APPEND objects
    hamilt_pw.cpp
    ks_pw/ekinetic_pw.cpp
    ks_pw/veff_pw.cpp
    ks_pw/nonlocal_pw.cpp
    ks_pw/meta_pw.cpp
    ks_pw/velocity_pw.cpp
)

if(ENABLE_LCAO)
    list(APPEND objects
        hamilt_lcao.cpp
        ks_lcao/operator_lcao.cpp
        ks_lcao/ekinetic_lcao.cpp
        ks_lcao/nonlocal_lcao.cpp
        ks_lcao/overlap_lcao.cpp
        ks_lcao/veff_lcao.cpp
        ks_lcao/meta_lcao.cpp
        ks_lcao/op_dftu_lcao.cpp
        ks_lcao/deepks_lcao.cpp
        ks_lcao/op_exx_lcao.cpp
    )
endif()

add_library(
    hamilt
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt)
endif()
