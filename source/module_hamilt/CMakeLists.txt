add_subdirectory(ks_pw)

add_library(
    hamilt
    OBJECT
    operator.cpp
    hamilt_pw.cpp
    hamilt_lcao.cpp
    ks_pw/ekinetic_pw.cpp
    ks_pw/veff_pw.cpp
    ks_pw/nonlocal_pw.cpp
    ks_pw/meta_pw.cpp
)
