add_library(
    xc_
    OBJECT
    xc_functional.cpp
    xc_functional_vxc.cpp
    xc_functional_gradcorr.cpp
    xc_functional_wrapper_xc.cpp
    xc_functional_wrapper_gcxc.cpp
    xc_functional_wrapper_tauxc.cpp
    xc_funct_exch_lda.cpp
    xc_funct_corr_lda.cpp
    xc_funct_exch_gga.cpp
    xc_funct_corr_gga.cpp
    xc_funct_hcth.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(xc_)
endif()

if(BUILD_TESTING)
if(ENABLE_LIBXC)
  add_subdirectory(test)
  add_subdirectory(kernels/test)
endif()
endif()
