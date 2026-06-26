remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
#include "module_/.h"
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_atom
  SOURCES sltk_atom_test.cpp ../sltk_atom.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_grid
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES sltk_grid_test.cpp ../sltk_grid.cpp ../sltk_atom.cpp ../atom_pack.cpp ../mpi_domain.cpp
    ../../../source_io/module_output/output.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_atom_arrange
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES sltk_atom_arrange_test.cpp  ../sltk_atom_arrange.cpp ../sltk_grid_driver.cpp ../sltk_grid.cpp
  ../sltk_atom.cpp ../atom_pack.cpp ../mpi_domain.cpp ../neighbor_rebuild_state.cpp
  ../../../source_io/module_output/output.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_mpi_domain
  LIBS ${math_libs}
  SOURCES mpi_domain_test.cpp ../mpi_domain.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_atom_pack
  LIBS parameter ${math_libs} base device cell_info
  SOURCES atom_pack_test.cpp ../atom_pack.cpp ../mpi_domain.cpp ../sltk_grid.cpp ../sltk_atom.cpp
    ../../../source_io/module_output/output.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_mpi_neighbor_proto
  LIBS ${math_libs}
  SOURCES mpi_neighbor_proto_test.cpp ../atom_pack.cpp ../mpi_domain.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_mpi_grid
  LIBS parameter ${math_libs} base device cell_info
  SOURCES mpi_grid_test.cpp ../sltk_grid_driver.cpp ../sltk_grid.cpp ../sltk_atom.cpp
    ../atom_pack.cpp ../mpi_domain.cpp ../neighbor_rebuild_state.cpp
    ../../../source_io/module_output/output.cpp
)

# Standalone performance probes used by the Phase 4 benchmark driver. They are
# regular build targets rather than CTest tests because benchmark repetitions
# and problem sizes are selected explicitly by the reporting script.
add_executable(MODULE_CELL_NEIGHBOR_benchmark
  neighbor_benchmark.cpp ../atom_pack.cpp ../mpi_domain.cpp)
target_link_libraries(MODULE_CELL_NEIGHBOR_benchmark
  parameter ${math_libs} base device Threads::Threads)
if(USE_OPENMP)
  target_link_libraries(MODULE_CELL_NEIGHBOR_benchmark OpenMP::OpenMP_CXX)
endif()

add_executable(MODULE_CELL_NEIGHBOR_mpi_benchmark
  mpi_neighbor_benchmark.cpp ../atom_pack.cpp ../mpi_domain.cpp)
target_link_libraries(MODULE_CELL_NEIGHBOR_mpi_benchmark
  parameter ${math_libs} base device Threads::Threads)
if(USE_OPENMP)
  target_link_libraries(MODULE_CELL_NEIGHBOR_mpi_benchmark OpenMP::OpenMP_CXX)
endif()

add_executable(MODULE_CELL_NEIGHBOR_grid_path_benchmark
  grid_path_benchmark.cpp ../atom_pack.cpp ../mpi_domain.cpp
  ../../../source_io/module_output/output.cpp)
target_link_libraries(MODULE_CELL_NEIGHBOR_grid_path_benchmark
  parameter ${math_libs} base device cell_info Threads::Threads)
if(USE_OPENMP)
  target_link_libraries(MODULE_CELL_NEIGHBOR_grid_path_benchmark OpenMP::OpenMP_CXX)
endif()
