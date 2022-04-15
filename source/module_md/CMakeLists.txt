add_library(
    md
    OBJECT
    MD_func.cpp
    run_md_classic.cpp
    LJ_potential.cpp
    DP_potential.cpp
    cmd_neighbor.cpp
	verlet.cpp
	NVE.cpp
    MSST.cpp
    NVT_ADS.cpp
    NVT_NHC.cpp
    FIRE.cpp
    Langevin.cpp
)

IF (BUILD_TESTING)
    add_subdirectory(test)
endif()