list(APPEND objects
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

add_library(
    xc
    OBJECT
    ${objects}
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()