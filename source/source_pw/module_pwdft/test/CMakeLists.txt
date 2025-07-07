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
	../../../source_cell/unitcell.cpp
	../../../source_io/output.cpp
	../../../source_cell/update_cell.cpp
	../../../source_cell/bcast_cell.cpp
	../../../source_cell/print_cell.cpp
	../../../source_cell/atom_spec.cpp
	../../../source_cell/atom_pseudo.cpp
	../../../source_cell/pseudo.cpp
	../../../source_cell/read_stru.cpp
	../../../source_cell/read_atom_species.cpp
	../../../source_cell/read_atoms.cpp
	../../../source_cell/read_pp.cpp
	../../../source_cell/read_pp_complete.cpp
	../../../source_cell/read_pp_upf100.cpp
	../../../source_cell/read_pp_upf201.cpp
	../../../source_cell/read_pp_vwr.cpp
	../../../source_cell/read_pp_blps.cpp
	../../../source_estate/read_pseudo.cpp
	../../../source_estate/cal_wfc.cpp
	../../../source_estate/cal_nelec_nband.cpp
	../../../source_estate/read_orb.cpp
)