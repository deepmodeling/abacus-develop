add_library(
    md
    OBJECT
    fire.cpp
    langevin.cpp
    md_base.cpp
    md_func.cpp
    msst.cpp
    nhchain.cpp
    run_md.cpp
    verlet.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(md)
endif()

if(BUILD_TESTING)
    add_subdirectory(test)
endif()
