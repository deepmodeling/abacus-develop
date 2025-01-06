remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

list(APPEND cell_source_files
    ../../../module_cell/update_cell.cpp
    ../../../module_cell/bcast_cell.cpp
    ../../../module_io/output.cpp
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
          ../../../module_io/orb_io.cpp
)

AddTest(
  TARGET bfgs_basic_test
  LIBS parameter ${math_libs} base device 
  SOURCES bfgs_basic_test.cpp ../bfgs_basic.cpp 
)


AddTest(
  TARGET bfgs_test
  LIBS parameter ${math_libs} base device 
  SOURCES bfgs_test.cpp ../bfgs.cpp ../ions_move_basic.cpp ${cell_source_files}
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
          ../../../module_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_cg_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_cg_test.cpp 
          ../ions_move_cg.cpp 
          ../ions_move_basic.cpp
          ../../../module_io/orb_io.cpp
          ${cell_source_files}
)

AddTest(
  TARGET ions_move_sd_test
  LIBS parameter ${math_libs} base device 
  SOURCES ions_move_sd_test.cpp ../ions_move_sd.cpp ../ions_move_basic.cpp ${cell_source_files}
)