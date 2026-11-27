# smoothgen is consumed by both NAO (atomic_radials) and PW
# (source_pw/module_proj/onsite_proj_init) code, so it must be built
# and linked regardless of ENABLE_LCAO.
add_library(nao_smoothgen OBJECT smoothgen.cpp)

if(ENABLE_COVERAGE)
  add_coverage(nao_smoothgen)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()

if(ENABLE_LCAO)
  add_library(
    numerical_atomic_orbitals OBJECT
    numerical_radial.cpp
    radial_set.cpp
    atomic_radials.cpp
    beta_radials.cpp
    sphbes_radials.cpp
    hydrogen_radials.cpp
    pswfc_radials.cpp
    radial_collection.cpp
    real_gaunt_table.cpp
    two_center_table.cpp
    two_center_integrator.cpp
    two_center_bundle.cpp)

  if(ENABLE_COVERAGE)
    add_coverage(numerical_atomic_orbitals)
  endif()

  if(BUILD_TESTING)
    if(ENABLE_MPI)
      add_subdirectory(test)
    endif()
  endif()
endif()
