add_library(
    orb
    OBJECT
    ORB_atomic.cpp
    ORB_atomic_lm.cpp
    ORB_control.cpp
    ORB_gaunt_table.cpp
    ORB_gen_tables.cpp
    ORB_nonlocal.cpp
    ORB_nonlocal_lm.cpp
    ORB_read.cpp
    ORB_table_alpha.cpp
    ORB_table_beta.cpp
    ORB_table_phi.cpp
    parallel_orbitals.cpp
)

set(CMAKE_CXX_STANDARD_REQUIRED ON) 
IF (BUILD_TESTING)
  set(CMAKE_CXX_STANDARD 14)
  add_subdirectory(test)
endif()