remove_definitions(-DUSE_PAW)
add_executable(
    test_deepks
    main_deepks.cpp klist_1.cpp LCAO_deepks_test_prep.cpp LCAO_deepks_test.cpp
    ../../../module_cell/unitcell.cpp
    ../../../module_cell/update_cell.cpp
    ../../../module_cell/bcast_cell.cpp
    ../../../module_cell/atom_spec.cpp
    ../../../module_cell/atom_pseudo.cpp
    ../../../module_cell/read_atoms.cpp
    ../../../module_cell/read_stru.cpp
    ../../../module_cell/read_atom_species.cpp
    ../../../module_cell/setup_nonlocal.cpp
    ../../../module_cell/pseudo.cpp
    ../../../module_cell/read_pp.cpp
    ../../../module_cell/read_pp_complete.cpp
    ../../../module_cell/read_pp_upf100.cpp
    ../../../module_cell/read_pp_upf201.cpp
    ../../../module_cell/read_pp_vwr.cpp
    ../../../module_cell/read_pp_blps.cpp
    ../../../module_hamilt_pw/hamilt_pwdft/soc.cpp
    ../../../module_io/output.cpp
    ../../../module_elecstate/read_pseudo.cpp
    ../../../module_elecstate/cal_nelec_nband.cpp
    ../../../module_elecstate/module_dm/density_matrix.cpp
    ../../../module_elecstate/module_dm/density_matrix_io.cpp
    ../../../module_hamilt_lcao/module_hcontainer/base_matrix.cpp 
    ../../../module_hamilt_lcao/module_hcontainer/hcontainer.cpp 
    ../../../module_hamilt_lcao/module_hcontainer/atom_pair.cpp
    ../../../module_hamilt_lcao/module_hcontainer/func_transfer.cpp
    ../../../module_hamilt_lcao/module_hcontainer/func_folding.cpp
    ../../../module_hamilt_lcao/module_hcontainer/transfer.cpp
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
