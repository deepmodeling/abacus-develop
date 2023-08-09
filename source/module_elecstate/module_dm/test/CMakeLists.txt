remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET read_dm_gamma
  LIBS ${math_libs} base device cell_info
  SOURCES test_readdmgamma.cpp ${ABACUS_SOURCE_DIR}/module_io/read_dm.cpp ${ABACUS_SOURCE_DIR}/module_io/write_dm.cpp ${ABACUS_SOURCE_DIR}/module_io/output.cpp
)

AddTest(
  TARGET io_dmk_io
  LIBS ${math_libs} base device cell_info symmetry
  SOURCES test_dmk_io.cpp ../dmk_io.cpp ${ABACUS_SOURCE_DIR}/module_io/output.cpp
  ${ABACUS_SOURCE_DIR}/module_cell/klist.cpp
  ${ABACUS_SOURCE_DIR}/module_cell/parallel_kpoints.cpp
)