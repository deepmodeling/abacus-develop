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
	../dfpt_pw_data.cpp
	../dfpt_kq_basis.cpp
	../../../source_cell/qlist.cpp
	../../../source_cell/reciprocal_grid.cpp
	../../../source_psi/psi.cpp
	# Plus_U test-support shim shared with the MPI-side dfpt tests.
	../test/dftu_test_support.cpp
)
