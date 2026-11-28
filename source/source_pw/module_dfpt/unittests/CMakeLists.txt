# Consolidated DFPT unit tests (the legacy test_serial/ and test/ directories
# were migrated here per the unittests/ layout rule; the module_stodft move in
# #8334 is the precedent). Target and CTest names are unchanged from the
# legacy layout, so the CI test matrix is unaffected.
#
# Serial targets drive the plane-wave FFT machinery directly (recip2real /
# real2recip on the shared rho/wfc grid). PW_Basis changes its member layout
# under __MPI, so the pw sources linked here must be compiled without __MPI
# in the same translation-unit ABI as the dfpt sources.
abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__MLALGO)

add_library(
  dfpt_planewave_serial OBJECT
  ../../../source_base/module_fft/fft_bundle.cpp
  ../../../source_base/module_fft/fft_cpu.cpp
  ../../../source_basis/module_pw/pw_basis.cpp
  ../../../source_basis/module_pw/pw_basis_k.cpp
  ../../../source_basis/module_pw/pw_basis_sup.cpp
  ../../../source_basis/module_pw/pw_distributeg.cpp
  ../../../source_basis/module_pw/pw_distg_method1.cpp
  ../../../source_basis/module_pw/pw_distg_method2.cpp
  ../../../source_basis/module_pw/pw_distributer.cpp
  ../../../source_basis/module_pw/pw_init.cpp
  ../../../source_basis/module_pw/pw_transform.cpp
  ../../../source_basis/module_pw/pw_transform_k.cpp)

addtest(
  TARGET MODULE_DFPT_pert_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES test_dfpt_pert.cpp
          ../dfpt_pert.cpp
          ../dfpt_pert_vkb.cpp
          ../dfpt_pert_nl.cpp
          ../dfpt_pw_data.cpp
          ../dfpt_pw_data_resp.cpp
          ../dfpt_kq_basis.cpp
          ../../../source_cell/qlist.cpp
          ../../../source_cell/reciprocal_grid.cpp
          ../../../source_psi/psi.cpp
          # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
          # (see test_dfpt_mocks.cpp).
          test_dfpt_mocks.cpp
          # Shared serial-side cell/basis/data fixture.
          test_dfpt_fixture.cpp)

addtest(
  TARGET MODULE_DFPT_rho_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES test_dfpt_rho.cpp
          test_dfpt_rho_mix.cpp
          ../dfpt_rho.cpp
          ../dfpt_rho_mix.cpp
          ../../../source_estate/module_charge/chg_tools.cpp
          ../dfpt_pw_data.cpp
          ../dfpt_pw_data_resp.cpp
          ../dfpt_kq_basis.cpp
          ../../../source_cell/qlist.cpp
          ../../../source_cell/reciprocal_grid.cpp
          ../../../source_psi/psi.cpp
          # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
          # (see test_dfpt_mocks.cpp).
          test_dfpt_mocks.cpp
          # Shared serial-side cell/basis/data fixture.
          test_dfpt_fixture.cpp)

addtest(
  TARGET MODULE_DFPT_phon_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES test_dfpt_phon_ewald.cpp
          test_dfpt_nscf.cpp
          test_dfpt_phon_elec.cpp
          test_dfpt_phon_asr.cpp
          test_dfpt_phon.cpp
          ../dfpt_phon.cpp
          ../dfpt_phon_asr.cpp
          ../dfpt_phon_ewald.cpp
          ../dfpt_phon_elec.cpp
          ../dfpt_nscf.cpp
          # Shared Ewald neighbor-shell generator (ion_ion real-space sum);
          # raw sources because hamilt_general is an MPI-ABI object library
          ../../../source_hamilt/module_ewald/h_ewald_pw.cpp
          ../../../source_hamilt/module_ewald/dnrm2.cpp
          ../dfpt_rho.cpp
          ../dfpt_rho_mix.cpp
          ../../../source_estate/module_charge/chg_tools.cpp
          ../dfpt_pert.cpp
          ../dfpt_pert_vkb.cpp
          ../dfpt_pert_nl.cpp
          ../dfpt_pw_data.cpp
          ../dfpt_pw_data_resp.cpp
          ../dfpt_kq_basis.cpp
          ../../../source_cell/qlist.cpp
          ../../../source_cell/reciprocal_grid.cpp
          ../../../source_psi/psi.cpp
          # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
          # (see test_dfpt_mocks.cpp).
          test_dfpt_mocks.cpp
          # Shared serial-side cell/basis/data fixture.
          test_dfpt_fixture.cpp)

addtest(
  TARGET MODULE_DFPT_q0_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES test_dfpt_q0.cpp
          test_dfpt_q0_pos.cpp
          ../dfpt_q0.cpp
          ../dfpt_q0_pos.cpp
          ../dfpt_pert.cpp
          ../dfpt_pert_vkb.cpp
          ../dfpt_pert_nl.cpp
          ../dfpt_pw_data.cpp
          ../dfpt_pw_data_resp.cpp
          ../dfpt_kq_basis.cpp
          ../../../source_cell/qlist.cpp
          ../../../source_cell/reciprocal_grid.cpp
          ../../../source_psi/psi.cpp
          # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
          # (see test_dfpt_mocks.cpp).
          test_dfpt_mocks.cpp
          # Shared serial-side cell/basis/data fixture.
          test_dfpt_fixture.cpp)

# QE-format dyn writer + the Gamma-first q-ordering validation consumed by
# the DFPT driver (test_dfpt_write_dyn.cpp matches dfpt_write_dyn.cpp,
# coding rule 11)
addtest(
  TARGET MODULE_DFPT_write_dyn_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES test_dfpt_write_dyn.cpp
          ../dfpt_write_dyn.cpp
          ../dfpt_pw_data.cpp
          ../dfpt_pw_data_resp.cpp
          ../../../source_cell/qlist.cpp
          ../../../source_cell/reciprocal_grid.cpp
          ../../../source_psi/psi.cpp
          # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
          # (see test_dfpt_mocks.cpp).
          test_dfpt_mocks.cpp
          # Shared serial-side cell/basis/data fixture.
          test_dfpt_fixture.cpp)

# Analytic libxc first-order XC kernel (empty TU without libxc, so the test
# only exists when the feature is compiled in)
if(ENABLE_LIBXC)
  addtest(
    TARGET MODULE_DFPT_xc_analytic_serial
    LIBS parameter dfpt_planewave_serial device base symmetry Libxc::xc
    SOURCES test_dfpt_xc_analytic.cpp
            ../dfpt_xc_analytic.cpp
            ../dfpt_rho.cpp
            ../dfpt_rho_mix.cpp
            ../../../source_estate/module_charge/chg_tools.cpp
            ../dfpt_pw_data.cpp
            ../dfpt_pw_data_resp.cpp
            ../dfpt_kq_basis.cpp
            ../../../source_cell/qlist.cpp
            ../../../source_cell/reciprocal_grid.cpp
            ../../../source_psi/psi.cpp
            # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
            # (see test_dfpt_mocks.cpp).
            test_dfpt_mocks.cpp
            # Shared serial-side cell/basis/data fixture.
            test_dfpt_fixture.cpp)
endif()

# MPI-build targets. The directory disables __MPI for the serial ABI above;
# these targets genuinely exercise the MPI feature, so they opt back in via
# KEEP_FEATURE_DEFINITIONS (see AddTest in cmake/Testing.cmake).
if(ENABLE_MPI)
  abacus_disable_feature_definitions(__MLALGO)
  abacus_disable_feature_definitions(__CUDA)
  abacus_disable_feature_definitions(__ROCM)
  abacus_disable_feature_definitions(__EXX)

  addtest(
    TARGET MODULE_DFPT_pw_data_test
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter base device symmetry
    SOURCES test_dfpt_pw_data.cpp
            ../dfpt_pw_data.cpp
            ../dfpt_pw_data_resp.cpp
            ../../../source_cell/qlist.cpp
            ../../../source_cell/reciprocal_grid.cpp
            ../../../source_psi/psi.cpp
            # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures.
            test_dfpt_mocks.cpp
            # Shared cubic-cell / stru_lib fixture.
            test_dfpt_stru_fixture.cpp)

  addtest(
    TARGET MODULE_DFPT_kq_basis_test
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter base device symmetry planewave
    SOURCES test_dfpt_kq_basis.cpp ../dfpt_kq_basis.cpp test_dfpt_mocks.cpp)

  addtest(
    TARGET MODULE_DFPT_stern_test
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter base device symmetry
    SOURCES test_dfpt_stern.cpp ../dfpt_stern.cpp test_dfpt_mocks.cpp)

  addtest(
    TARGET MODULE_DFPT_pw_run_test
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter base device symmetry planewave
    SOURCES test_dfpt_pw_run.cpp
            ../dfpt_pw.cpp
            ../dfpt_pw_init.cpp
            ../dfpt_pw_run.cpp
            ../dfpt_pw_solve.cpp
            ../dfpt_pw_q0.cpp
            ../dfpt_pw_data.cpp
            ../dfpt_pw_data_resp.cpp
            ../dfpt_pert.cpp
            ../dfpt_pert_vkb.cpp
            ../dfpt_pert_nl.cpp
            ../dfpt_kq_basis.cpp
            ../dfpt_stern.cpp
            ../dfpt_nscf.cpp
            ../dfpt_rho.cpp
            ../dfpt_rho_mix.cpp
            ../../../source_estate/module_charge/chg_tools.cpp
            ../dfpt_phon.cpp
            ../dfpt_phon_asr.cpp
            ../dfpt_phon_ewald.cpp
            ../dfpt_phon_elec.cpp
            # Shared Ewald neighbor-shell generator (ion_ion real-space sum);
            # raw sources because hamilt_general is an MPI-ABI object library
            ../../../source_hamilt/module_ewald/h_ewald_pw.cpp
            ../../../source_hamilt/module_ewald/dnrm2.cpp
            ../dfpt_q0.cpp
            ../dfpt_q0_pos.cpp
            ../dfpt_metal.cpp
            ../dfpt_hamilt_shift.cpp
            ../../../source_cell/qlist.cpp
            ../../../source_cell/reciprocal_grid.cpp
            ../../../source_psi/psi.cpp
            # Shared ctor/dtor stubs for the cell/spepot/charge/Plus_U closures.
            test_dfpt_mocks.cpp
            # Shared cubic-cell / stru_lib fixture.
            test_dfpt_stru_fixture.cpp)
endif()

# The little-group enumeration is self-contained (no cell fixture), so its
# test links only the module itself plus base.
addtest(
  TARGET MODULE_DFPT_symm_serial
  LIBS device base
  SOURCES test_dfpt_symm.cpp ../dfpt_symm.cpp)
