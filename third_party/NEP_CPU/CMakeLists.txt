cmake_minimum_required(VERSION 3.16)
project(NEP_CPU CXX)

set(NEP_CPU_SOURCES
  src/nep.cpp
  src/ewald_nep.cpp
  src/neighbor_nep.cpp
)

add_library(nep SHARED ${NEP_CPU_SOURCES})
target_include_directories(nep PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(nep PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  POSITION_INDEPENDENT_CODE ON
)

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(nep PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS nep
  EXPORT NEPTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

install(FILES
  src/nep.h
  src/ewald_nep.h
  src/neighbor_nep.h
  src/nep_utilities.h
  src/dftd3para.h
  DESTINATION include
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/NEPConfigVersion.cmake"
  VERSION 1.0.0
  COMPATIBILITY AnyNewerVersion
)

install(EXPORT NEPTargets
  FILE NEPTargets.cmake
  NAMESPACE NEP::
  DESTINATION lib/cmake/NEP
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/NEPConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/NEPConfig.cmake"
  @ONLY
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/NEPConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/NEPConfigVersion.cmake"
  DESTINATION lib/cmake/NEP
)
