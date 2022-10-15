list(APPEND objects
    hamilt_pw.cpp
    ks_pw/ekinetic_pw.cpp
    ks_pw/veff_pw.cpp
    ks_pw/nonlocal_pw.cpp
    ks_pw/meta_pw.cpp
)

if(ENABLE_LCAO)
    list(APPEND objects
        hamilt_lcao.cpp
    )
endif()

add_library(
    hamilt
    OBJECT
    ${objects}
)

add_subdirectory(ks_pw)
