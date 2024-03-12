remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

AddTest(
  TARGET lattice_change_methods_test
  LIBS ${math_libs} base device
  SOURCES lattice_change_methods_test.cpp ../lattice_change_methods.cpp ../lattice_change_basic.cpp
)

AddTest(
  TARGET lattice_change_basic_test
  LIBS ${math_libs} base device
  SOURCES lattice_change_basic_test.cpp ../lattice_change_basic.cpp
)

AddTest(
  TARGET lattice_change_cg_test
  LIBS ${math_libs} base device
  SOURCES lattice_change_cg_test.cpp ../lattice_change_cg.cpp ../lattice_change_basic.cpp
)

AddTest(
  TARGET bfgs_basic_test
  LIBS ${math_libs} base device
  SOURCES bfgs_basic_test.cpp ../bfgs_basic.cpp
)

AddTest(
  TARGET ions_move_methods_test
  LIBS ${math_libs} base device
  SOURCES ions_move_methods_test.cpp ../ions_move_methods.cpp ../ions_move_basic.cpp ../ions_move_bfgs.cpp ../ions_move_cg.cpp ../ions_move_sd.cpp ../bfgs_basic.cpp
)

AddTest(
  TARGET ions_move_basic_test
  LIBS ${math_libs} base device
  SOURCES ions_move_basic_test.cpp ../ions_move_basic.cpp
)

AddTest(
  TARGET ions_move_bfgs_test
  LIBS ${math_libs} base device
  SOURCES ions_move_bfgs_test.cpp ../ions_move_bfgs.cpp ../ions_move_basic.cpp ../bfgs_basic.cpp
)

AddTest(
  TARGET ions_move_cg_test
  LIBS ${math_libs} base device
  SOURCES ions_move_cg_test.cpp ../ions_move_cg.cpp ../ions_move_basic.cpp
)

AddTest(
  TARGET ions_move_sd_test
  LIBS ${math_libs} base device
  SOURCES ions_move_sd_test.cpp ../ions_move_sd.cpp ../ions_move_basic.cpp
)