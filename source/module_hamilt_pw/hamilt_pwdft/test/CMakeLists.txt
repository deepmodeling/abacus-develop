remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET pwdft_soc
  LIBS parameter  ${math_libs}
  SOURCES soc_test.cpp ../soc.cpp
  	../../../source_base/global_variable.cpp
  	../../../source_base/global_function.cpp
	../../../source_base/global_file.cpp
	../../../source_base/tool_quit.cpp
	../../../source_base/global_variable.cpp
	../../../source_base/global_file.cpp
	../../../source_base/memory.cpp
	../../../source_base/timer.cpp
	../../../source_base/blas_connector_base.cpp ../../../source_base/blas_connector_vector.cpp ../../../source_base/blas_connector_matrix.cpp
	../../../source_base/parallel_global.cpp
	../../../source_base/parallel_comm.cpp
	../../../source_base/parallel_common.cpp
	../../../source_base/parallel_reduce.cpp
)

AddTest(
	TARGET radial_proj_test
	LIBS parameter  base device ${math_libs}
	SOURCES radial_proj_test.cpp ../radial_proj.cpp
)

AddTest(
	TARGET structure_factor_test
	LIBS parameter ${math_libs} base device planewave 
	SOURCES structure_factor_test.cpp ../structure_factor.cpp ../parallel_grid.cpp
	../../../module_cell/unitcell.cpp
	../../../module_io/output.cpp
	../../../module_cell/update_cell.cpp
	../../../module_cell/bcast_cell.cpp
	../../../module_cell/print_cell.cpp
	../../../module_cell/atom_spec.cpp
	../../../module_cell/atom_pseudo.cpp
	../../../module_cell/pseudo.cpp
	../../../module_cell/read_stru.cpp
	../../../module_cell/read_atom_species.cpp
	../../../module_cell/read_atoms.cpp
	../../../module_cell/read_pp.cpp
	../../../module_cell/read_pp_complete.cpp
	../../../module_cell/read_pp_upf100.cpp
	../../../module_cell/read_pp_upf201.cpp
	../../../module_cell/read_pp_vwr.cpp
	../../../module_cell/read_pp_blps.cpp
	../../../module_elecstate/read_pseudo.cpp
	../../../module_elecstate/cal_wfc.cpp
	../../../module_elecstate/cal_nelec_nband.cpp
	../../../module_elecstate/read_orb.cpp
)