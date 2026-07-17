abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
#include "module_/.h"
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_atom
  SOURCES sltk_atom_test.cpp ../sltk_atom.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_grid
  LIBS parameter base device cell_info
  SOURCES sltk_grid_test.cpp ../sltk_grid.cpp ../sltk_atom.cpp

)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_sltk_atom_arrange
  LIBS parameter base device cell_info
  SOURCES sltk_atom_arrange_test.cpp  ../sltk_atom_arrange.cpp ../sltk_grid_driver.cpp ../sltk_grid.cpp
  ../sltk_atom.cpp

)