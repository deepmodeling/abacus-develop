add_library(
    md
    OBJECT
    MD_func.cpp
    run_md.cpp
	mdrun.cpp
	verlet.cpp
    MSST.cpp
    Nose_Hoover.cpp
    FIRE.cpp
    Langevin.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(md)
endif()

if(BUILD_TESTING)
    add_subdirectory(test)
endif()
