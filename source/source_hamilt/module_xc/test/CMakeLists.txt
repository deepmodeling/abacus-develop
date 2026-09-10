AddTest(
  TARGET MODULE_HAMILT_XCTest_PBE
  LIBS parameter MPI::MPI_CXX Libxc::xc
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_lda_wrap.cpp ../xc_gga_wrap.cpp ../xc_gga_corr.cpp ../xc_lda_corr.cpp ../xc_gga_exch.cpp ../xc_lda_exch.cpp ../xc_hcth.cpp
  ../libxc_gga_wrap.cpp ../libxc_setup.cpp
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_HSE
  LIBS parameter MPI::MPI_CXX Libxc::xc
  SOURCES test_xc1.cpp ../xc_functional.cpp ../libxc_setup.cpp
)


AddTest(
  TARGET MODULE_HAMILT_XCTest_PZ_SPN
  LIBS parameter MPI::MPI_CXX Libxc::xc
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_lda_wrap.cpp ../xc_gga_wrap.cpp ../xc_gga_corr.cpp ../xc_lda_corr.cpp ../xc_gga_exch.cpp ../xc_lda_exch.cpp ../xc_hcth.cpp
  ../libxc_gga_wrap.cpp ../libxc_lda_wrap.cpp ../libxc_setup.cpp
)

if (USE_CUDA)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_cuda.cpp)
endif()
if (USE_ROCM)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_rocm.cpp)
endif()
AddTest(
  TARGET MODULE_HAMILT_XCTest_GRADCORR
  LIBS parameter MPI::MPI_CXX Libxc::xc psi device container
  SOURCES test_xc3.cpp ../xc_grad.cpp ../xc_grad_prepare.cpp ../xc_grad_kernel.cpp ../xc_grad_assemble.cpp ../xc_grad_wfc.cpp ../xc_grad_utils.cpp ../xc_functional.cpp
  ../xc_lda_wrap.cpp ../xc_gga_wrap.cpp
  ../libxc_setup.cpp
  ../libxc_lda_wrap.cpp
  ../libxc_gga_wrap.cpp
  ../libxc_mgga_wrap.cpp
  ../xc_gga_corr.cpp ../xc_lda_corr.cpp ../xc_gga_exch.cpp
  ../xc_lda_exch.cpp ../xc_hcth.cpp
  ../../../source_base/matrix.cpp
  ../../../source_base/memory_recorder.cpp
  ../../../source_base/global_variable.cpp
  ../../../source_base/libm/branred.cpp
  ../../../source_base/libm/sincos.cpp
  ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
  ../../../source_base/module_fft/fft_bundle.cpp
  ../../../source_base/module_fft/fft_cpu.cpp
  ${FFT_SRC}
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_SCAN
  LIBS parameter MPI::MPI_CXX Libxc::xc
  SOURCES test_xc4.cpp ../xc_functional.cpp ../xc_lda_wrap.cpp
    ../xc_gga_wrap.cpp
    ../libxc_setup.cpp
    ../libxc_lda_wrap.cpp
    ../libxc_gga_wrap.cpp
    ../libxc_mgga_wrap.cpp
    ../xc_gga_corr.cpp ../xc_lda_corr.cpp
    ../xc_gga_exch.cpp ../xc_lda_exch.cpp ../xc_hcth.cpp
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_VXC
  LIBS parameter MPI::MPI_CXX Libxc::xc psi device container
  SOURCES test_xc5.cpp test_libxc_tools.cpp ../xc_grad.cpp ../xc_grad_prepare.cpp ../xc_grad_kernel.cpp ../xc_grad_assemble.cpp ../xc_grad_wfc.cpp ../xc_grad_utils.cpp ../xc_functional.cpp
    ../xc_lda_wrap.cpp ../xc_gga_wrap.cpp
    ../libxc_setup.cpp
    ../libxc_lda_wrap.cpp
    ../libxc_gga_wrap.cpp
    ../libxc_mgga_wrap.cpp
    ../xc_gga_corr.cpp ../xc_lda_corr.cpp ../xc_gga_exch.cpp
    ../xc_lda_exch.cpp ../xc_hcth.cpp
    ../xc_pot.cpp
    ../xc_functional_ncgga_sf.cpp
    ../xc_ncgga_radial.cpp
    ../libxc_pot.cpp
    ../libxc_tools.cpp
    ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
    ../../../source_base/matrix.cpp
    ../../../source_base/memory_recorder.cpp
  ../../../source_base/global_variable.cpp
    ../../../source_base/timer.cpp
    ../../../source_base/libm/branred.cpp
    ../../../source_base/libm/sincos.cpp
    ../../../source_base/module_fft/fft_bundle.cpp
    ../../../source_base/module_fft/fft_cpu.cpp
    ${FFT_SRC}
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_NCGGA_RADIAL
  SOURCES test_xc_ncgga_radial.cpp ../xc_ncgga_radial.cpp
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD
  LIBS xc_ planewave parameter MPI::MPI_CXX Libxc::xc psi device container base
  SOURCES test_xc_functional_ncgga_sf.cpp
)
set_tests_properties(
  MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1"
)
add_test(
  NAME MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD_MPI2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD>
          ${MPIEXEC_POSTFLAGS}
)
set_tests_properties(
  MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD_MPI2
  PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1"
    PROCESSORS 2
    WORKING_DIRECTORY $<TARGET_FILE_DIR:MODULE_HAMILT_XCTest_NCGGA_DISCRETE_FD>
)

AddTest(
  TARGET MODULE_HAMILT_XCTest_SCANL_LAPL
  LIBS parameter MPI::MPI_CXX Libxc::xc
  SOURCES test_xc6.cpp ../xc_functional.cpp ../xc_lda_wrap.cpp
    ../xc_gga_wrap.cpp
    ../libxc_setup.cpp
    ../libxc_lda_wrap.cpp
    ../libxc_gga_wrap.cpp
    ../libxc_mgga_wrap.cpp
    ../xc_gga_corr.cpp ../xc_lda_corr.cpp
    ../xc_gga_exch.cpp ../xc_lda_exch.cpp ../xc_hcth.cpp
)

if (USE_CUDA)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_cuda.cpp)
endif()
if (USE_ROCM)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_rocm.cpp)
endif()
AddTest(
  TARGET MODULE_HAMILT_XCTest_LAPL
  LIBS parameter MPI::MPI_CXX Libxc::xc psi device container
  SOURCES test_xc7.cpp ../xc_grad.cpp ../xc_grad_prepare.cpp ../xc_grad_kernel.cpp ../xc_grad_assemble.cpp ../xc_grad_wfc.cpp ../xc_grad_utils.cpp ../xc_functional.cpp
  ../xc_lda_wrap.cpp ../xc_gga_wrap.cpp
  ../libxc_setup.cpp
  ../libxc_lda_wrap.cpp
  ../libxc_gga_wrap.cpp
  ../libxc_mgga_wrap.cpp
  ../xc_gga_corr.cpp ../xc_lda_corr.cpp ../xc_gga_exch.cpp
  ../xc_lda_exch.cpp ../xc_hcth.cpp
  ../../../source_base/matrix.cpp
  ../../../source_base/memory_recorder.cpp
  ../../../source_base/global_variable.cpp
  ../../../source_base/libm/branred.cpp
  ../../../source_base/libm/sincos.cpp
  ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
  ../../../source_base/module_fft/fft_bundle.cpp
  ../../../source_base/module_fft/fft_cpu.cpp
  ${FFT_SRC}
)

# gradcorr dispatches the noncollinear stress path even in focused FFT tests.
foreach(xc_gradient_test MODULE_HAMILT_XCTest_GRADCORR MODULE_HAMILT_XCTest_LAPL)
  target_sources(${xc_gradient_test} PRIVATE
    ../xc_functional_ncgga_sf.cpp ../xc_ncgga_radial.cpp
    ../libxc_pot.cpp ../libxc_tools.cpp ../../../source_base/timer.cpp)
endforeach()
