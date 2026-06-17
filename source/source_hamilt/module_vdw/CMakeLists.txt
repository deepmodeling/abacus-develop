set(vdw_sources
    vdwd2_parameters.cpp
    vdwd3_parameters_tab.cpp
    vdwd3_parameters.cpp
    vdwd2.cpp
    vdwd3.cpp
    vdwd3_autoset_xcname.cpp
    vdwd3_autoset_xcparam.cpp
    vdw.cpp
)

if(ENABLE_DFTD4)
  list(APPEND vdw_sources vdwd4.cpp)
endif()

add_library(vdw OBJECT ${vdw_sources})

if(ENABLE_DFTD4)
  target_link_libraries(vdw PUBLIC dftd4::dftd4)
  target_compile_definitions(vdw PUBLIC __DFTD4)
endif()

if(ENABLE_COVERAGE)
  add_coverage(vdw)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
