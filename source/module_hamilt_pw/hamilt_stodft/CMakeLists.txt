list(APPEND hamilt_stodft_srcs
    sto_iter.cpp
    sto_hchi.cpp
    sto_wf.cpp
    sto_func.cpp
    sto_forces.cpp
    sto_stress_pw.cpp
)

add_library(
    hamilt_stodft
    OBJECT
    ${hamilt_stodft_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_stodft)
endif()
