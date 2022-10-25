add_library(
    operator_ks_pw
    OBJECT
    ekinetic_pw.cpp
    veff_pw.cpp
    nonlocal_pw.cpp
    meta_pw.cpp
    velocity_pw.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_pw)
endif()
