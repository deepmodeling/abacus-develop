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
    module_charge/charge.cpp
    module_charge/charge_broyden.cpp
    module_charge/charge_extra.cpp
    module_charge/charge_mixing.cpp
    module_charge/charge_pulay.cpp
    module_charge/symmetry_rho.cpp
    energy.cpp
    magnetism.cpp
    occupy.cpp
)

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
