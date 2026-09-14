abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_DFPT_pw_data_test
  LIBS parameter base device symmetry
  SOURCES dfpt_pw_data_test.cpp
  	../dfpt_pw_data.cpp
  	../../../source_cell/qlist.cpp
  	../../../source_cell/reciprocal_grid.cpp
  	../../../source_psi/psi.cpp
  	# Shared ctor/dtor stubs for the cell/spepot/Plus_U link closures.
  	dfpt_test_mocks.cpp
  	# Shared cubic-cell / stru_lib fixture.
  	dfpt_stru_fixture.cpp
)

AddTest(
  TARGET MODULE_DFPT_kq_basis_test
  LIBS parameter base device symmetry planewave
  SOURCES dfpt_kq_basis_test.cpp
    ../dfpt_kq_basis.cpp
    dfpt_test_mocks.cpp
)

AddTest(
  TARGET MODULE_DFPT_stern_test
  LIBS parameter base device symmetry
  SOURCES dfpt_stern_test.cpp
  	../dfpt_stern.cpp
  	dfpt_test_mocks.cpp
)

AddTest(
  TARGET MODULE_DFPT_pw_run_test
  LIBS parameter base device symmetry planewave
  SOURCES dfpt_pw_run_test.cpp
	../dfpt_pw.cpp
	../dfpt_pw_init.cpp
	../dfpt_pw_run.cpp
	../dfpt_pw_solve.cpp
	../dfpt_pw_q0.cpp
	../dfpt_pw_data.cpp
	../dfpt_pert.cpp
	../dfpt_pert_vkb.cpp
	../dfpt_pert_nl.cpp
	../dfpt_kq_basis.cpp
	../dfpt_stern.cpp
	../dfpt_rho.cpp
	../dfpt_phon.cpp
	../dfpt_phon_ewald.cpp
	../dfpt_phon_elec.cpp
	../dfpt_q0.cpp
	../dfpt_q0_pos.cpp
	../dfpt_metal.cpp
	../dfpt_hamilt_shift.cpp
	../../../source_cell/qlist.cpp
	../../../source_cell/reciprocal_grid.cpp
	../../../source_psi/psi.cpp
	# Shared ctor/dtor stubs for the cell/spepot/charge/Plus_U closures.
	dfpt_test_mocks.cpp
	# Shared cubic-cell / stru_lib fixture.
	dfpt_stru_fixture.cpp
)
