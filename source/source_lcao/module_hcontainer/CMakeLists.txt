if(ENABLE_LCAO)

list(APPEND objects
    base_matrix.cpp
    atom_pair.cpp
    hcontainer.cpp
    output_hcontainer.cpp
    func_folding.cpp
    transfer.cpp
    func_transfer.cpp
)

add_library(
    hcontainer
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hcontainer)
endif()

IF (BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()

endif()