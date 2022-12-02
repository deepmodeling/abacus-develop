list(APPEND objects
    elecstate.cpp
    elecstate_pw.cpp
    elecstate_pw_sdft.cpp
    potentials/gatefield.cpp
    potentials/efield.cpp
    potentials/H_Hartree_pw.cpp
    potentials/pot_xc.cpp
    potentials/pot_local.cpp
    potentials/potential_new.cpp
    potentials/write_pot.cpp
    src/elecstate_multi_device.cpp
)

if (USE_CUDA) 
  list(APPEND objects src/cuda/elecstate_multi_device.cu)
elseif(USE_ROCM)
  list(APPEND objects src/rocm/elecstate_multi_device.cu)
endif()

if(ENABLE_LCAO)
  list(APPEND objects
      elecstate_lcao.cpp
      elecstate_lcao_tddft.cpp
      potentials/H_TDDFT_pw.cpp
  )
endif()

add_library(
    elecstate
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(elecstate)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
