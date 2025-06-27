remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)
remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET relax_new_line_search
  LIBS parameter 
  SOURCES line_search_test.cpp ../line_search.cpp ../../source_base/global_variable.cpp ../../source_base/global_file.cpp ../../source_base/global_function.cpp ../../source_base/memory.cpp ../../source_base/timer.cpp ../../source_base/tool_quit.cpp
)

AddTest(
  TARGET relax_new_relax
  SOURCES  relax_test.cpp ../relax_sync.cpp ../line_search.cpp ../../source_base/tool_quit.cpp ../../source_base/global_variable.cpp ../../source_base/global_file.cpp ../../source_base/memory.cpp ../../source_base/timer.cpp
    ../../source_base/matrix3.cpp ../../source_base/intarray.cpp ../../source_base/tool_title.cpp
    ../../source_base/global_function.cpp ../../source_base/complexmatrix.cpp ../../source_base/matrix.cpp
    ../../source_base/complexarray.cpp ../../source_base/tool_quit.cpp ../../source_base/realarray.cpp
    ../../source_base/blas_connector_base.cpp ../../source_base/blas_connector_vector.cpp ../../source_base/blas_connector_matrix.cpp
    ../../source_cell/update_cell.cpp ../../source_cell/print_cell.cpp  ../../source_cell/bcast_cell.cpp ../../source_io/output.cpp
    LIBS parameter ${math_libs} 
)

list(APPEND cell_source_files
    ../../source_cell/update_cell.cpp
    ../../source_cell/bcast_cell.cpp
    ../../source_cell/print_cell.cpp
    ../../source_io/output.cpp
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
          ../../source_io/orb_io.cpp
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
          ../../source_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_cg_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_cg_test.cpp 
          ../ions_move_cg.cpp 
          ../ions_move_basic.cpp
          ../../source_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_sd_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_sd_test.cpp ../ions_move_sd.cpp ../ions_move_basic.cpp ${cell_source_files}
)