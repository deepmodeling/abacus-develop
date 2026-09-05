set(vdw_sources
    vdw_xcname.cpp
    vdwd2_parameters.cpp
    vdwd2.cpp
    vdwd3_data.cpp
    vdwd3_parameters.cpp
    vdwd3_evaluator.cpp
    vdwd3.cpp
    vdw.cpp
)

if(ENABLE_DFTD4)
  list(APPEND vdw_sources vdwd4.cpp)
endif()

add_library(vdw OBJECT ${vdw_sources})

if(ENABLE_DFTD4)
  target_link_libraries(vdw PUBLIC dftd4::dftd4)
endif()

if(ENABLE_COVERAGE)
  add_coverage(vdw)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
