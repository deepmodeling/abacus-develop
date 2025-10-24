list(APPEND objects
    elecstate.cpp
    elecstate_energy_terms.cpp
    elecstate_energy.cpp
    elecstate_exx.cpp
    elecstate_print.cpp
    elecstate_tools.cpp
    elecstate_pw.cpp
    elecstate_pw_sdft.cpp
    elecstate_pw_cal_tau.cpp
    module_pot/gatefield.cpp
    module_pot/efield.cpp
    module_pot/H_Hartree_pw.cpp
    module_pot/pot_xc.cpp
    module_pot/pot_local.cpp
    module_pot/pot_local_paw.cpp
    module_pot/potential_new.cpp
    module_pot/potential_types.cpp
    module_pot/pot_sep.cpp
    module_charge/chgmixing.cpp
    module_charge/charge.cpp
    module_charge/charge_init.cpp
    module_charge/charge_mpi.cpp
    module_charge/charge_extra.cpp
    module_charge/charge_mixing.cpp
    module_charge/charge_mixing_dmr.cpp
    module_charge/charge_mixing_residual.cpp
    module_charge/charge_mixing_preconditioner.cpp
    module_charge/charge_mixing_rho.cpp
    module_charge/charge_mixing_uspp.cpp
    module_charge/symmetry_rho.cpp
    module_charge/symmetry_rhog.cpp
    fp_energy.cpp
    magnetism.cpp
    occupy.cpp
    cal_ux.cpp
    read_orb.cpp
    cal_nelec_nband.cpp
    read_pseudo.cpp
    cal_wfc.cpp
    setup_estate_pw.cpp
    update_pot.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      elecstate_lcao.cpp
      elecstate_lcao_cal_tau.cpp
      module_pot/H_TDDFT_pw.cpp
      module_dm/setup_dm.cpp
      module_dm/density_matrix.cpp
      module_dm/density_matrix_io.cpp
      module_dm/cal_dm_psi.cpp
      module_dm/cal_edm_tddft.cpp
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
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(test_mpi)
  endif()
endif()

if(ENABLE_LCAO)
  add_subdirectory(module_dm)
endif()
