remove_definitions(-D__MPI -D__LCAO )
add_definitions(-D__NORMAL)

list(APPEND depend_files 
  ../md_func.cpp
  ../../source_cell/unitcell.cpp
  ../../source_cell/update_cell.cpp
  ../../source_cell/bcast_cell.cpp
  ../../source_cell/print_cell.cpp
  ../../source_cell/atom_spec.cpp
  ../../source_cell/read_stru.cpp
  ../../source_cell/read_atom_species.cpp
  ../../source_cell/atom_pseudo.cpp
  ../../source_cell/read_atoms.cpp
  ../../source_cell/pseudo.cpp
  ../../source_cell/read_pp.cpp
  ../../source_cell/read_pp_complete.cpp
  ../../source_cell/read_pp_upf100.cpp
  ../../source_cell/read_pp_upf201.cpp
  ../../source_cell/read_pp_vwr.cpp
  ../../source_cell/read_pp_blps.cpp
  ../../source_base/matrix3.cpp
  ../../source_base/matrix.cpp
  ../../source_base/timer.cpp
  ../../source_base/module_external/blas_connector_base.cpp
  ../../source_base/module_external/blas_connector_matrix.cpp
  ../../source_base/module_external/blas_connector_vector.cpp
  ../../source_base/memory.cpp
  ../../source_base/global_variable.cpp
  ../../source_base/global_function.cpp
  ../../source_base/global_file.cpp
  ../../source_base/tool_title.cpp
  ../../source_base/tool_check.cpp
  ../../source_base/tool_quit.cpp
  ../../source_base/intarray.cpp
  ../../source_base/realarray.cpp
  ../../source_base/complexarray.cpp
  ../../source_base/complexmatrix.cpp
  ../../source_base/global_variable.cpp
  ../../source_base/libm/branred.cpp
  ../../source_base/libm/sincos.cpp
  ../../source_base/math_integral.cpp
  ../../source_cell/module_neighbor/sltk_atom_arrange.cpp
  ../../source_cell/module_neighbor/sltk_atom.cpp
  ../../source_cell/module_neighbor/sltk_grid.cpp
  ../../source_cell/module_neighbor/sltk_grid_driver.cpp
  ../../source_io/output.cpp
  ../../source_io/output_log.cpp
  ../../source_io/print_info.cpp
  ../../source_io/cif_io.cpp
  ../../source_esolver/esolver_lj.cpp
  ../../source_base/parallel_reduce.cpp
  ../../source_base/parallel_global.cpp
  ../../source_base/parallel_comm.cpp
  ../../source_estate/read_pseudo.cpp
  ../../source_estate/cal_wfc.cpp
  ../../source_estate/cal_nelec_nband.cpp
  ../../source_estate/read_orb.cpp
)

AddTest(
  TARGET MODULE_MD_LJ_pot
  LIBS parameter ${math_libs} psi device 
  SOURCES lj_pot_test.cpp  
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_func
  LIBS parameter ${math_libs} psi device 
  SOURCES md_func_test.cpp  
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_fire
  LIBS parameter ${math_libs} psi device 
  SOURCES fire_test.cpp
  ../md_base.cpp
  ../fire.cpp
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_verlet
  LIBS parameter ${math_libs} psi device 
  SOURCES verlet_test.cpp
  ../md_base.cpp
  ../verlet.cpp
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_nhc
  LIBS parameter ${math_libs} psi device 
  SOURCES nhchain_test.cpp
  ../md_base.cpp
  ../nhchain.cpp
  ../../source_io/output.cpp
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_msst
  LIBS parameter ${math_libs} psi device 
  SOURCES msst_test.cpp
  ../md_base.cpp
  ../msst.cpp
  ../../source_cell/update_cell.cpp
  ../../source_io/output.cpp
  ${depend_files}
)

AddTest(
  TARGET MODULE_MD_lgv
  LIBS parameter ${math_libs} psi device 
  SOURCES langevin_test.cpp
  ../md_base.cpp
  ../langevin.cpp
  ${depend_files}
)
