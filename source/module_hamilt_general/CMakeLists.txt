add_subdirectory(module_vdw)
add_subdirectory(module_surchem)
add_subdirectory(module_xc)

list(APPEND objects
    operator.cpp
    module_ewald/H_Ewald_pw.cpp
    module_ewald/dnrm2.cpp
)

add_library(
    hamilt_general
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_general)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
