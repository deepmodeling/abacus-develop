list(APPEND objects
      dftu.cpp
      dftu_force.cpp
      dftu_yukawa.cpp
      dftu_folding.cpp
      dftu_io.cpp
      dftu_tools.cpp
      dftu_occup.cpp
      dftu_hamilt.cpp
      dftu_pw.cpp
)

add_library(
    dftu
    OBJECT
    ${objects}
)

# Workaround for NVC++ internal compiler error in dftu_yukawa.cpp
# Use -O1 optimization level for this specific file
if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    set_source_files_properties(dftu_yukawa.cpp PROPERTIES COMPILE_FLAGS "-O1")
endif()

if(ENABLE_COVERAGE)
  add_coverage(dftu)
endif()
