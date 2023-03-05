add_library(
        vdw
        OBJECT
        vdwd2_parameters.cpp
        vdwd3_parameters_tab.cpp
        vdwd3_parameters.cpp
        vdwd2.cpp
        vdwd3.cpp
        vdw.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(vdw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()