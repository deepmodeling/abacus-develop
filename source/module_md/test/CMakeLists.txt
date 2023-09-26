remove_definitions(-D__MPI -D__LCAO )
add_definitions(-D__NORMAL)
remove_definitions(-DUSE_PAW)

list(APPEND depend_files 
  ../md_func.cpp
  ../../module_io/input.cpp
  ../../module_cell/unitcell.cpp
  ../../module_cell/atom_spec.cpp
  ../../module_cell/atom_pseudo.cpp
  ../../module_cell/read_atoms.cpp
  ../../module_cell/read_cell_pseudopots.cpp
  ../../module_cell/pseudo.cpp
  ../../module_cell/read_pp.cpp
  ../../module_cell/read_pp_upf100.cpp
  ../../module_cell/read_pp_upf201.cpp
  ../../module_cell/read_pp_vwr.cpp
  ../../module_cell/read_pp_blps.cpp
  ../../module_base/matrix3.cpp
  ../../module_base/matrix.cpp
  ../../module_base/timer.cpp
  ../../module_base/memory.cpp
  ../../module_base/global_variable.cpp
  ../../module_base/global_function.cpp
  ../../module_base/global_file.cpp
  ../../module_base/tool_title.cpp
  ../../module_base/tool_check.cpp
  ../../module_base/tool_quit.cpp
  ../../module_base/intarray.cpp
  ../../module_base/realarray.cpp
  ../../module_base/complexarray.cpp
  ../../module_base/complexmatrix.cpp
  ../../module_base/global_variable.cpp
  ../../module_base/libm/branred.cpp
  ../../module_base/libm/sincos.cpp
  ../../module_cell/module_neighbor/sltk_adjacent_set.cpp
  ../../module_cell/module_neighbor/sltk_atom_arrange.cpp
  ../../module_cell/module_neighbor/sltk_atom.cpp
  ../../module_cell/module_neighbor/sltk_atom_input.cpp
  ../../module_cell/module_neighbor/sltk_grid.cpp
  ../../module_cell/module_neighbor/sltk_grid_driver.cpp
  ../../module_io/output.cpp
  ../../module_io/print_info.cpp
  ../../module_esolver/esolver_lj.cpp
  ../../module_base/parallel_reduce.cpp
  ../../module_base/parallel_global.cpp
)

AddTest(
  TARGET md_LJ_pot
  LIBS ${math_libs} psi device formatter
  SOURCES lj_pot_test.cpp  
  ${depend_files}
)

AddTest(
  TARGET md_func
  LIBS ${math_libs} psi device formatter
  SOURCES md_func_test.cpp  
  ${depend_files}
)

AddTest(
  TARGET md_fire
  LIBS ${math_libs} psi device formatter
  SOURCES fire_test.cpp
  ../md_base.cpp
  ../fire.cpp
  ${depend_files}
)

AddTest(
  TARGET md_verlet
  LIBS ${math_libs} psi device formatter
  SOURCES verlet_test.cpp
  ../md_base.cpp
  ../verlet.cpp
  ${depend_files}
)

AddTest(
  TARGET md_nhc
  LIBS ${math_libs} psi device formatter
  SOURCES nhchain_test.cpp
  ../md_base.cpp
  ../nhchain.cpp
  ${depend_files}
)

AddTest(
  TARGET md_msst
  LIBS ${math_libs} psi device formatter
  SOURCES msst_test.cpp
  ../md_base.cpp
  ../msst.cpp
  ${depend_files}
)

AddTest(
  TARGET md_lgv
  LIBS ${math_libs} psi device formatter
  SOURCES langevin_test.cpp
  ../md_base.cpp
  ../langevin.cpp
  ${depend_files}
)
