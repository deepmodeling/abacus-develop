list(APPEND objects
    elecstate.cpp
    elecstate_getters.cpp
    elecstate_energy_terms.cpp
    elecstate_energy.cpp
    elecstate_exx.cpp
    elecstate_print.cpp
    elecstate_pw.cpp
    elecstate_pw_sdft.cpp
    potentials/gatefield.cpp
    potentials/efield.cpp
    potentials/H_Hartree_pw.cpp
    potentials/pot_xc.cpp
    potentials/pot_local.cpp
    potentials/pot_local_paw.cpp
    potentials/potential_new.cpp
    potentials/potential_types.cpp
    module_charge/charge.cpp
    module_charge/charge_init.cpp
    module_charge/charge_mpi.cpp
    module_charge/charge_extra.cpp
    module_charge/charge_mixing.cpp
    module_charge/symmetry_rho.cpp
    module_charge/symmetry_rhog.cpp
    fp_energy.cpp
    magnetism.cpp
    occupy.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      elecstate_lcao.cpp
      elecstate_lcao_tddft.cpp
      potentials/H_TDDFT_pw.cpp
      module_dm/density_matrix.cpp
      module_dm/cal_dm_psi.cpp
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
  if(ENABLE_MPI)
    add_subdirectory(test_mpi)
  endif()
endif()

if(ENABLE_LCAO)
  add_subdirectory(module_dm)
endif()