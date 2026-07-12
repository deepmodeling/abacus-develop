abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_RELAX_relax_new_line_search
  LIBS parameter
  SOURCES line_search_test.cpp ../line_search.cpp ../../source_base/global_variable.cpp ../../source_base/global_file.cpp ../../source_base/global_function.cpp ../../source_base/memory_recorder.cpp ../../source_base/timer.cpp ../../source_base/tool_quit.cpp
)

AddTest(
  TARGET MODULE_RELAX_relax_new_relax
  SOURCES  relax_test.cpp ../relax_sync.cpp ../line_search.cpp ../../source_base/tool_quit.cpp ../../source_base/global_variable.cpp ../../source_base/global_file.cpp ../../source_base/memory_recorder.cpp ../../source_base/timer.cpp
    ../../source_base/matrix3.cpp ../../source_base/intarray.cpp ../../source_base/tool_title.cpp
    ../../source_base/global_function.cpp ../../source_base/complexmatrix.cpp ../../source_base/matrix.cpp
    ../../source_base/complexarray.cpp ../../source_base/tool_quit.cpp ../../source_base/realarray.cpp
    ../../source_base/module_external/blas_connector_base.cpp ../../source_base/module_external/blas_connector_vector.cpp ../../source_base/module_external/blas_connector_matrix.cpp
    ../../source_cell/update_cell.cpp ../../source_cell/print_cell.cpp  ../../source_cell/bcast_cell.cpp ../../source_base/output.cpp
    LIBS parameter
)

list(APPEND cell_source_files
    ../../source_cell/update_cell.cpp
    ../../source_cell/bcast_cell.cpp
    ../../source_cell/print_cell.cpp
)
AddTest(
  TARGET MODULE_RELAX_lattice_change_methods_test
  LIBS parameter base device
  SOURCES lattice_change_methods_test.cpp ../lattice_change_methods.cpp ../lattice_change_basic.cpp ../relax_data.cpp mock_remake_cell.cpp
)

AddTest(
  TARGET MODULE_RELAX_lattice_change_basic_test
  LIBS parameter base device
  SOURCES lattice_change_basic_test.cpp ../lattice_change_basic.cpp ../relax_data.cpp mock_remake_cell.cpp
)

AddTest(
  TARGET MODULE_RELAX_lattice_change_cg_test
  LIBS parameter base device
  SOURCES lattice_change_cg_test.cpp
          ../lattice_change_cg.cpp
          ../lattice_change_basic.cpp
          ../cg_base.cpp
          ../relax_data.cpp
          mock_remake_cell.cpp
          ../../source_io/module_output/orb_io.cpp
)

AddTest(
  TARGET MODULE_RELAX_bfgs_basic_test
  LIBS parameter base device symmetry
  SOURCES bfgs_basic_test.cpp ../bfgs_basic.cpp ../relax_data.cpp ../ions_move_basic.cpp
          ../../source_io/module_output/orb_io.cpp
          ${cell_source_files}
)


AddTest(
  TARGET MODULE_RELAX_bfgs_test
  LIBS parameter base device
  SOURCES bfgs_test.cpp ../ions_move_bfgs2.cpp ../ions_move_basic.cpp ../matrix_methods.cpp ../relax_data.cpp ${cell_source_files}
)

AddTest(
  TARGET MODULE_RELAX_ions_move_basic_test
  LIBS parameter base device
  SOURCES ions_move_basic_test.cpp ../ions_move_basic.cpp ../relax_data.cpp ${cell_source_files}
)

AddTest(
  TARGET MODULE_RELAX_ions_move_bfgs_test
  LIBS parameter base device
  SOURCES ions_move_bfgs_test.cpp 
          ../ions_move_bfgs.cpp 
          ../ions_move_basic.cpp 
          ../bfgs_basic.cpp
          ../relax_data.cpp
          ../../source_io/module_output/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET MODULE_RELAX_ions_move_methods_test
  LIBS parameter base device
  SOURCES ions_move_methods_test.cpp
          ../ions_move_methods.cpp
          ../ions_move_bfgs.cpp
          ../ions_move_bfgs2.cpp
          ../ions_move_cg.cpp
          ../ions_move_sd.cpp
          ../ions_move_lbfgs.cpp
          ../ions_move_basic.cpp
          ../bfgs_basic.cpp
          ../cg_base.cpp
          ../matrix_methods.cpp
          ../relax_data.cpp
          ../../source_io/module_output/orb_io.cpp
          ../../source_cell/bcast_cell.cpp
          ../../source_cell/print_cell.cpp
)

AddTest(
  TARGET MODULE_RELAX_ions_move_cg_test
  LIBS parameter base device
  SOURCES ions_move_cg_test.cpp 
          ../ions_move_cg.cpp 
          ../cg_base.cpp
          ../ions_move_basic.cpp
          ../relax_data.cpp
          ../../source_io/module_output/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET MODULE_RELAX_ions_move_sd_test
  LIBS parameter base device
  SOURCES ions_move_sd_test.cpp ../ions_move_sd.cpp ../ions_move_basic.cpp ../relax_data.cpp ${cell_source_files}
)
