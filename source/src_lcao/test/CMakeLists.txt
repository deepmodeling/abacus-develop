remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)
AddTest(
  TARGET hsolver_gamma_rho
  LIBS ${math_libs} base neighbor 
  SOURCES gamma_rho_test.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/ft.cpp ../../src_pw/pw_basis.cpp ../../src_pw/pw_complement.cpp ../../src_pw/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../gint_gamma.cpp ../gint_gamma_rho.cpp ../gint_tools.cpp ../grid_technique.cpp ../grid_base_beta.cpp ../grid_meshk.cpp ../grid_meshcell.cpp ../grid_bigcell.cpp ../grid_meshball.cpp ../../module_orbital/ORB_atomic.cpp ../../module_orbital/ORB_atomic_lm.cpp ../../module_orbital/ORB_gaunt_table.cpp ../../module_orbital/ORB_table_phi.cpp ../../module_orbital/ORB_table_beta.cpp ../../module_orbital/ORB_table_alpha.cpp ../../module_orbital/ORB_nonlocal_lm.cpp ../../module_orbital/ORB_nonlocal.cpp ../../module_cell/atom_spec.cpp ../../module_cell/atom_pseudo.cpp ../../module_cell/pseudo_nc.cpp ../../src_io/output.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
