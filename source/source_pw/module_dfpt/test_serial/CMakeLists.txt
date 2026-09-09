abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__MLALGO)

# The DFPT first-order potential kernels (dfpt_pert.cpp) drive the plane-wave
# FFT machinery directly (recip2real / real2recip on the shared rho/wfc grid).
# PW_Basis changes its member layout under __MPI, so the pw sources linked
# here must be compiled without __MPI in the same translation-unit ABI as the
# dfpt sources (mirrors module_pw/test_serial).
add_library(
    dfpt_planewave_serial
    OBJECT
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
    ../../../source_basis/module_pw/pw_transform_k.cpp
)

AddTest(
  TARGET MODULE_DFPT_pert_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES dfpt_pert_serial_test.cpp
	../dfpt_pert.cpp
	../dfpt_pert_vkb.cpp
	../dfpt_pert_nl.cpp
	../dfpt_pw_data.cpp
	../dfpt_kq_basis.cpp
	../../../source_cell/qlist.cpp
	../../../source_cell/reciprocal_grid.cpp
	../../../source_psi/psi.cpp
	# Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
  # (see test/dfpt_test_mocks.cpp).
	../test/dfpt_test_mocks.cpp
	# Shared serial-side cell/basis/data fixture.
	dfpt_serial_fixture.cpp
)

AddTest(
  TARGET MODULE_DFPT_rho_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES dfpt_rho_serial_test.cpp
  ../dfpt_rho.cpp
  ../dfpt_pw_data.cpp
  ../dfpt_kq_basis.cpp
  ../../../source_cell/qlist.cpp
  ../../../source_cell/reciprocal_grid.cpp
  ../../../source_psi/psi.cpp
  # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
  # (see test/dfpt_test_mocks.cpp).
  ../test/dfpt_test_mocks.cpp
  # Shared serial-side cell/basis/data fixture.
  dfpt_serial_fixture.cpp
)

AddTest(
  TARGET MODULE_DFPT_phon_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES dfpt_phon_serial_test.cpp
  ../dfpt_phon.cpp
  ../dfpt_phon_ewald.cpp
  ../dfpt_phon_elec.cpp
  ../dfpt_pert.cpp
  ../dfpt_pert_vkb.cpp
  ../dfpt_pert_nl.cpp
  ../dfpt_pw_data.cpp
  ../dfpt_kq_basis.cpp
  ../../../source_cell/qlist.cpp
  ../../../source_cell/reciprocal_grid.cpp
  ../../../source_psi/psi.cpp
  # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
  # (see test/dfpt_test_mocks.cpp).
  ../test/dfpt_test_mocks.cpp
  # Shared serial-side cell/basis/data fixture.
  dfpt_serial_fixture.cpp
)

AddTest(
  TARGET MODULE_DFPT_q0_serial
  LIBS parameter dfpt_planewave_serial device base symmetry
  SOURCES dfpt_q0_serial_test.cpp
  ../dfpt_q0.cpp
  ../dfpt_q0_pos.cpp
  ../dfpt_pert.cpp
  ../dfpt_pert_vkb.cpp
  ../dfpt_pert_nl.cpp
  ../dfpt_pw_data.cpp
  ../dfpt_kq_basis.cpp
  ../../../source_cell/qlist.cpp
  ../../../source_cell/reciprocal_grid.cpp
  ../../../source_psi/psi.cpp
  # Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures
  # (see test/dfpt_test_mocks.cpp).
  ../test/dfpt_test_mocks.cpp
  # Shared serial-side cell/basis/data fixture.
  dfpt_serial_fixture.cpp
)
