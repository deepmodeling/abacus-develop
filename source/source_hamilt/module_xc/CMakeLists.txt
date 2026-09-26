add_library(
    xc_
    OBJECT
    xc_functional.cpp
    xc_pot.cpp
    xc_grad.cpp
    xc_grad_prepare.cpp
    xc_grad_kernel.cpp
    xc_grad_assemble.cpp
    xc_grad_wfc.cpp
    xc_grad_utils.cpp
    xc_rvv10.cpp
    xc_rvv10_pw.cpp
    xc_lda_wrap.cpp
    xc_gga_wrap.cpp
    xc_lda_exch.cpp
    xc_lda_corr.cpp
    xc_gga_exch.cpp
    xc_gga_corr.cpp
    xc_hcth.cpp
    libxc_setup.cpp
    libxc_tools.cpp
    libxc_pot.cpp
    libxc_lda_wrap.cpp
    libxc_gga_wrap.cpp
    libxc_mgga_wrap.cpp
    exx_info.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(xc_)
endif()

if(BUILD_TESTING)
  if(NOT ENABLE_MPI AND NOT USE_CUDA AND NOT USE_ROCM)
    AddTest(
      TARGET MODULE_HAMILT_XCTest_RVV10_PW
      LIBS xc_ planewave base parameter device
      SOURCES test/test_xc_rvv10_pw.cpp
    )
    if(ENABLE_LIBXC)
      target_link_libraries(MODULE_HAMILT_XCTest_RVV10_PW PRIVATE Libxc::xc)
    endif()
    set_tests_properties(MODULE_HAMILT_XCTest_RVV10_PW PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  endif()
  AddTest(
    TARGET MODULE_HAMILT_XCTest_RVV10_CORE
    SOURCES test/test_xc_rvv10.cpp xc_rvv10.cpp ../../source_base/cubic_spline.cpp
  )
if(ENABLE_MPI)
if(ENABLE_LIBXC)
  add_subdirectory(test)
  add_subdirectory(kernels/test)
endif()
endif()
endif()
