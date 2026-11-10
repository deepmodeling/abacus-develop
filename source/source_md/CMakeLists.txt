add_library(
  md OBJECT
  fire.cpp
  langevin.cpp
  md_base.cpp
  md_func.cpp
  msst.cpp
  nhchain.cpp
  run_md.cpp
  verlet.cpp)

add_subdirectory(module_snap)
add_subdirectory(module_lj)

if(ENABLE_COVERAGE)
  add_coverage(md)
endif()

add_subdirectory(module_plumed)

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(unittests)
  endif()
endif()
