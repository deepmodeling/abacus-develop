remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)
remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET relax_new_line_search
  LIBS parameter 
  SOURCES line_search_test.cpp ../line_search.cpp ../../module_base/global_variable.cpp ../../module_base/global_file.cpp ../../module_base/global_function.cpp ../../module_base/memory.cpp ../../module_base/timer.cpp ../../module_base/tool_quit.cpp
)

AddTest(
  TARGET relax_new_relax
  SOURCES  relax_test.cpp ../relax_sync.cpp ../line_search.cpp ../../module_base/tool_quit.cpp ../../module_base/global_variable.cpp ../../module_base/global_file.cpp ../../module_base/memory.cpp ../../module_base/timer.cpp
    ../../module_base/matrix3.cpp ../../module_base/intarray.cpp ../../module_base/tool_title.cpp
    ../../module_base/global_function.cpp ../../module_base/complexmatrix.cpp ../../module_base/matrix.cpp
    ../../module_base/complexarray.cpp ../../module_base/tool_quit.cpp ../../module_base/realarray.cpp
    ../../module_base/blas_connector_base.cpp ../../module_base/blas_connector_vector.cpp ../../module_base/blas_connector_matrix.cpp
    ../../module_cell/update_cell.cpp ../../module_cell/print_cell.cpp  ../../module_cell/bcast_cell.cpp ../../module_io/output.cpp
    LIBS parameter ${math_libs} 
)

list(APPEND cell_source_files
    ../../module_cell/update_cell.cpp
    ../../module_cell/bcast_cell.cpp
    ../../module_cell/print_cell.cpp
    ../../module_io/output.cpp
)
AddTest(
  TARGET lattice_change_methods_test
  LIBS parameter ${math_libs} base device 
  SOURCES lattice_change_methods_test.cpp ../lattice_change_methods.cpp ../lattice_change_basic.cpp
)

AddTest(
  TARGET lattice_change_basic_test
  LIBS parameter ${math_libs} base device 
  SOURCES lattice_change_basic_test.cpp ../lattice_change_basic.cpp
)

AddTest(
  TARGET lattice_change_cg_test
  LIBS parameter ${math_libs} base device 
  SOURCES lattice_change_cg_test.cpp 
          ../lattice_change_cg.cpp 
          ../lattice_change_basic.cpp
          ../../module_io/orb_io.cpp
)

AddTest(
  TARGET bfgs_basic_test
  LIBS parameter ${math_libs} base device 
  SOURCES bfgs_basic_test.cpp ../bfgs_basic.cpp 
)


AddTest(
  TARGET bfgs_test
  LIBS parameter ${math_libs} base device 
  SOURCES bfgs_test.cpp ../bfgs.cpp ../ions_move_basic.cpp ../matrix_methods.cpp ${cell_source_files}
)

AddTest(
  TARGET ions_move_basic_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_basic_test.cpp ../ions_move_basic.cpp ${cell_source_files}
)

AddTest(
  TARGET ions_move_bfgs_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_bfgs_test.cpp 
          ../ions_move_bfgs.cpp 
          ../ions_move_basic.cpp 
          ../bfgs_basic.cpp
          ../../module_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_cg_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_cg_test.cpp 
          ../ions_move_cg.cpp 
          ../ions_move_basic.cpp
          ../../module_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_sd_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_sd_test.cpp ../ions_move_sd.cpp ../ions_move_basic.cpp ${cell_source_files}
)