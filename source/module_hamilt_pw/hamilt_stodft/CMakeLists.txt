list(APPEND hamilt_stodft_srcs
    hamilt_sdft_pw.cpp
    sto_iter.cpp
    sto_che.cpp
    sto_wf.cpp
    sto_func.cpp
    sto_forces.cpp
    sto_stress_pw.cpp
    sto_tool.cpp
    sto_elecond.cpp
    sto_dos.cpp
)

add_library(
    hamilt_stodft
    OBJECT
    ${hamilt_stodft_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_stodft)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()