remove_definitions(-DUSE_PAW)

AddTest(
  TARGET XCTest_PBE
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBEsol
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_BP
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_revPBE
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_WC
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_BLYP
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_OLYP
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_HCTH
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PWLDA
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PZ
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_SLATER1
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_SLATER_RXC
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PW
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_LYP
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_VWN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_WIGNER
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_HL
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_GL
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PW91
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_OPTX
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_WCX
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBE0
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBE_LibXC
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_HSE
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_SCAN0
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_KSDT
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_KT2
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_R2SCAN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_LB07
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_BMK
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_HF
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc1.cpp ../xc_functional.cpp
)

AddTest(
  TARGET XCTest_PZ_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBE_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_BP_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_revPBE_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_SLATER1_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_SLATER_RXC_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_P86_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBE0_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBEsol_SPN
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PBE_SPN_LibXC
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_PZ_SPN_LibXC
  LIBS MPI::MPI_CXX Libxc::xc # required by global.h; for details, `remove_definitions(-D__MPI)`.
  SOURCES test_xc2.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
)

AddTest(
  TARGET XCTest_GRADCORR
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container
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
  TARGET XCTest_GRADWFC
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container
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
  LIBS MPI::MPI_CXX Libxc::xc
  SOURCES test_xc4.cpp ../xc_functional.cpp ../xc_functional_wrapper_xc.cpp
    ../xc_functional_wrapper_gcxc.cpp ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp 
    ../xc_funct_exch_gga.cpp ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp 
)

AddTest(
  TARGET XCTest_VXC
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container formatter
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

AddTest(
  TARGET XCTest_VXC_Libxc
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container formatter
  SOURCES test_xc5.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../../../module_base/matrix.cpp ../../../module_base/timer.cpp
    ../xc_functional_vxc.cpp
    ../../../module_base/matrix.cpp
    ../../../module_base/memory.cpp
    ../../../module_base/timer.cpp
    ../../../module_base/libm/branred.cpp
    ../../../module_base/libm/sincos.cpp
)

AddTest(
  TARGET XCTest_VXC_meta
  LIBS MPI::MPI_CXX Libxc::xc ${math_libs} psi device container formatter
  SOURCES test_xc5.cpp ../xc_functional_gradcorr.cpp ../xc_functional.cpp
    ../xc_functional_wrapper_xc.cpp ../xc_functional_wrapper_gcxc.cpp
    ../xc_functional_wrapper_tauxc.cpp
    ../xc_funct_corr_gga.cpp ../xc_funct_corr_lda.cpp ../xc_funct_exch_gga.cpp
    ../xc_funct_exch_lda.cpp ../xc_funct_hcth.cpp
    ../../../module_base/matrix.cpp ../../../module_base/timer.cpp
    ../xc_functional_vxc.cpp
    ../../../module_base/matrix.cpp
    ../../../module_base/memory.cpp
    ../../../module_base/timer.cpp
    ../../../module_base/libm/branred.cpp
    ../../../module_base/libm/sincos.cpp
)

