remove_definitions(-DUSE_PAW)

AddTest(
  TARGET XCTest_PBE
  LIBS MPI::MPI_CXX Libxc::xc parameter # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_HSE
  LIBS MPI::MPI_CXX Libxc::xc parameter # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)


AddTest(
  TARGET XCTest_PZ_SPN
  LIBS MPI::MPI_CXX Libxc::xc parameter # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)


AddTest(
  TARGET XCTest_GRADCORR
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container parameter
  SOURCES test_xc3.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../../../module_base/matrix.cpp
    ../../../module_base/memory.cpp
    ../../../module_base/libm/branred.cpp
    ../../../module_base/libm/sincos.cpp
)

AddTest(
  TARGET XCTest_SCAN
  LIBS MPI::MPI_CXX Libxc::xc parameter
  SOURCES test_xc4.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp
    ../xc_functional_wrapper_gcxc.cpp ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp 
    ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp 
)

AddTest(
  TARGET XCTest_VXC
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container parameter
  SOURCES test_xc5.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../xc_functional_vxc.cpp
    ../../../module_base/matrix.cpp
    ../../../module_base/memory.cpp
    ../../../module_base/timer.cpp
    ../../../module_base/libm/branred.cpp
    ../../../module_base/libm/sincos.cpp
)
