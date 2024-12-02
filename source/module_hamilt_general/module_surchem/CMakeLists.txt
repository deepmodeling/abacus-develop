add_library(
    surchem
    OBJECT
    surchem.cpp
    H_correction_pw.cpp
    cal_epsilon.cpp
    cal_pseudo.cpp
    cal_totn.cpp
    cal_vcav.cpp
    cal_vel.cpp
    minimize_cg.cpp
    sol_force.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(surchem)
endif()

IF (BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
