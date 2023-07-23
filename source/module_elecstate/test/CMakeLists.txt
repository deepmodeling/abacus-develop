remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__UT_USE_CUDA)
remove_definitions(-D__UT_USE_ROCM)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
remove_definitions(-D_OPENMP)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET Elecstate_Op_UTs
  LIBS ${math_libs} psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
)

AddTest(
  TARGET elecstate_occupy
  LIBS ${math_libs} base device
  SOURCES elecstate_occupy_test.cpp ../occupy.cpp
)

AddTest(
  TARGET elecstate_magnetism
  LIBS ${math_libs} base device
  SOURCES elecstate_magnetism_test.cpp ../magnetism.cpp
)

AddTest(
  TARGET elecstate_fp_energy
  LIBS ${math_libs} base device
  SOURCES elecstate_fp_energy_test.cpp ../fp_energy.cpp
)

AddTest(
  TARGET elecstate_print
  LIBS ${math_libs} base device
  SOURCES elecstate_print_test.cpp ../elecstate_print.cpp ../occupy.cpp
)

AddTest(
  TARGET elecstate_base
  LIBS ${math_libs} base device
  SOURCES elecstate_base_test.cpp ../elecstate.cpp ../occupy.cpp ../../module_psi/psi.cpp
)

AddTest(
  TARGET elecstate_pw
  LIBS ${math_libs} planewave_serial base device
  SOURCES elecstate_pw_test.cpp
    ../elecstate_pw.cpp
    ../elecstate.cpp
    ../occupy.cpp
    ../../module_psi/psi.cpp
    ../../module_psi/kernels/memory_op.cpp
)

AddTest(
  TARGET elecstate_energy
  LIBS ${math_libs} base device
  SOURCES elecstate_energy_test.cpp ../elecstate_energy.cpp ../fp_energy.cpp
)

AddTest(
  TARGET potentials_base
  SOURCES potentials_base_test.cpp
)

AddTest(
  TARGET potentials_new
  LIBS ${math_libs} base device
  SOURCES potential_new_test.cpp ../potentials/potential_new.cpp
)

AddTest(
  TARGET charge_test
  LIBS ${math_libs} planewave_serial base device cell_simple
  SOURCES charge_test.cpp ../module_charge/charge.cpp
    ../../module_io/output.cpp
)

list(APPEND cell_simple_srcs
    ../../module_cell/unitcell.cpp
    ../../module_cell/read_atoms.cpp
    ../../module_cell/read_cell_pseudopots.cpp
    ../../module_cell/atom_spec.cpp
    ../../module_cell/atom_pseudo.cpp
    ../../module_cell/pseudo_nc.cpp
    ../../module_cell/read_pp.cpp
    ../../module_cell/read_pp_upf201.cpp
    ../../module_cell/read_pp_upf100.cpp
    ../../module_cell/read_pp_vwr.cpp
    ../../module_cell/read_pp_blps.cpp
)

AddTest(
  TARGET charge_mixing
  LIBS  base ${math_libs}  psi  device planewave_serial
  SOURCES charge_mixing_test.cpp  
  ../module_charge/charge_mixing.cpp
  ../module_charge/charge_broyden.cpp
  ../module_charge/charge_pulay.cpp
  ../../module_io/output.cpp
  ${cell_simple_srcs}
)

AddTest(
  TARGET charge_extra
  LIBS ${math_libs} base device
  SOURCES charge_extra_test.cpp ../module_charge/charge_extra.cpp ../../module_io/read_cube.cpp ../../module_io/write_cube.cpp
  ../../module_io/output.cpp ${cell_simple_srcs}
)