add_library(dftb_native OBJECT skf_data.cpp sk_block.cpp sk_matrix.cpp eigensolver.cpp charge_mixing.cpp periodic_scc.cpp)
target_link_libraries(dftb_native PRIVATE abacus::compile_requirements)

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(test)
endif()
