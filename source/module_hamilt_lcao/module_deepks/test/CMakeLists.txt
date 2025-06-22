add_executable(
    test_deepks
    main_deepks.cpp klist_1.cpp LCAO_deepks_test_prep.cpp LCAO_deepks_test.cpp
    ../../../source_cell/unitcell.cpp
    ../../../source_cell/update_cell.cpp
    ../../../source_cell/bcast_cell.cpp
    ../../../source_cell/atom_spec.cpp
    ../../../source_cell/atom_pseudo.cpp
    ../../../source_cell/read_atoms.cpp
    ../../../source_cell/read_stru.cpp
    ../../../source_cell/print_cell.cpp
    ../../../source_cell/read_atom_species.cpp
    ../../../source_cell/setup_nonlocal.cpp
    ../../../source_cell/pseudo.cpp
    ../../../source_cell/read_pp.cpp
    ../../../source_cell/read_pp_complete.cpp
    ../../../source_cell/read_pp_upf100.cpp
    ../../../source_cell/read_pp_upf201.cpp
    ../../../source_cell/read_pp_vwr.cpp
    ../../../source_cell/read_pp_blps.cpp
    ../../../module_hamilt_pw/hamilt_pwdft/soc.cpp
    ../../../module_io/output.cpp
    ../../../module_io/sparse_matrix.cpp
    ../../../module_elecstate/read_pseudo.cpp
    ../../../module_elecstate/cal_wfc.cpp
    ../../../module_elecstate/read_orb.cpp
    ../../../module_elecstate/cal_nelec_nband.cpp
    ../../../module_elecstate/module_dm/density_matrix.cpp
    ../../../module_elecstate/module_dm/density_matrix_io.cpp
    ../../../module_hamilt_lcao/module_hcontainer/base_matrix.cpp 
    ../../../module_hamilt_lcao/module_hcontainer/hcontainer.cpp 
    ../../../module_hamilt_lcao/module_hcontainer/atom_pair.cpp
    ../../../module_hamilt_lcao/module_hcontainer/func_transfer.cpp
    ../../../module_hamilt_lcao/module_hcontainer/func_folding.cpp
    ../../../module_hamilt_lcao/module_hcontainer/transfer.cpp
    ../../../module_hamilt_lcao/module_hcontainer/output_hcontainer.cpp
    ../../../module_hamilt_lcao/hamilt_lcaodft/operator_lcao/deepks_lcao.cpp
    ../../../module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.cpp
    ../../../module_hamilt_general/operator.cpp
)

target_link_libraries(
    test_deepks
    base device parameter deepks psi planewave neighbor container
    orb gint numerical_atomic_orbitals
    ${math_libs}
)

if(ENABLE_COVERAGE)
  add_coverage(test_deepks)
endif()


install(
    TARGETS test_deepks
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../../../tests
)
