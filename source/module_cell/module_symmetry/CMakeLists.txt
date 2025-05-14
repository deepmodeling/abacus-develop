add_library(
    symmetry
    OBJECT
    symm_other.cpp
    symmetry_basic.cpp
    symm_analysis.cpp
    symm_check.cpp
    symm_getgroup.cpp
    symm_hermite.cpp
    symm_lattice.cpp
    symm_magnetic.cpp
    symm_pricell.cpp
    symm_rho.cpp
    symmetry.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(symmetry)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
