AddTest(
  TARGET XCTest_PBE
  LIBS parameter MPI::MPI_CXX Libxc::xc  # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
  ../xc_functional_libxc_wrapper_gcxc.cpp ../xc_functional_libxc.cpp
)

AddTest(
  TARGET XCTest_HSE
  LIBS parameter MPI::MPI_CXX Libxc::xc  # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp ../xc_functional_libxc.cpp
)


AddTest(
  TARGET XCTest_PZ_SPN
  LIBS parameter MPI::MPI_CXX Libxc::xc  # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
  ../xc_functional_libxc_wrapper_gcxc.cpp ../xc_functional_libxc_wrapper_xc.cpp ../xc_functional_libxc.cpp
)

if (USE_CUDA)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_cuda.cpp)
endif()
if (USE_ROCM)
list(APPEND FFT_SRC ../../../source_base/module_fft/fft_rocm.cpp)
endif()
AddTest(
  TARGET XCTest_GRADCORR
  LIBS parameter MPI::MPI_CXX Libxc::xc ${math_libs} psi device container 
  SOURCES test_xc3.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_libxc.cpp
    ../xc_functional_libxc_wrapper_xc.cpp
    ../xc_functional_libxc_wrapper_gcxc.cpp
    ../xc_functional_libxc_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../../../source_base/matrix.cpp
    ../../../source_base/memory.cpp
    ../../../source_base/libm/branred.cpp
    ../../../source_base/libm/sincos.cpp
    ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
    ../../../source_base/module_fft/fft_bundle.cpp
    ../../../source_base/module_fft/fft_cpu.cpp
    ${FFT_SRC}
)

AddTest(
  TARGET XCTest_SCAN
  LIBS parameter MPI::MPI_CXX Libxc::xc 
  SOURCES test_xc4.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp
    ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_libxc.cpp
    ../xc_functional_libxc_wrapper_xc.cpp
    ../xc_functional_libxc_wrapper_gcxc.cpp
    ../xc_functional_libxc_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp 
    ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp 
)

AddTest(
  TARGET XCTest_VXC
  LIBS parameter MPI::MPI_CXX Libxc::xc ${math_libs} psi device container 
  SOURCES test_xc5.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_libxc.cpp
    ../xc_functional_libxc_wrapper_xc.cpp
    ../xc_functional_libxc_wrapper_gcxc.cpp
    ../xc_functional_libxc_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../xc_functional_vxc.cpp
    ../xc_functional_libxc_vxc.cpp
    ../xc_functional_libxc_tools.cpp
    ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
    ../../../source_base/matrix.cpp
    ../../../source_base/memory.cpp
    ../../../source_base/timer.cpp
    ../../../source_base/libm/branred.cpp
    ../../../source_base/libm/sincos.cpp
    ../../../source_base/module_fft/fft_bundle.cpp
    ../../../source_base/module_fft/fft_cpu.cpp
    ${FFT_SRC}
)