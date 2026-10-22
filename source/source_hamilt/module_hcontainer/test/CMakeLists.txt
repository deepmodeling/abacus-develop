if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_LCAO_hcontainer_test
    LIBS parameter psi base device
    SOURCES test_hcontainer.cpp ../base_matrix.cpp ../hcontainer.cpp
            ../atom_pair.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp)

  addtest(
    TARGET MODULE_LCAO_hcontainer_complex_test
    LIBS parameter psi base device
    SOURCES test_hcontainer_complex.cpp ../base_matrix.cpp ../hcontainer.cpp
            ../atom_pair.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp)

  # This driver only measures time and memory; keep it with the opt-in benchmarks.
  if(ENABLE_GOOGLEBENCH)
    addtest(
      TARGET MODULE_LCAO_hcontainer_cost_test
      LIBS parameter psi base device
      SOURCES test_hcontainer_time.cpp
              ../base_matrix.cpp
              ../hcontainer.cpp
              ../atom_pair.cpp
              ../../../source_basis/module_ao/parallel_orbitals.cpp
              tmp_mocks.cpp)
  endif()

  addtest(
    TARGET MODULE_LCAO_hcontainer_folding_test
    LIBS parameter psi base device
    SOURCES test_func_folding.cpp
            ../base_matrix.cpp
            ../hcontainer.cpp
            ../atom_pair.cpp
            ../func_folding.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp
            tmp_mocks.cpp)

  addtest(
    TARGET MODULE_LCAO_hcontainer_transfer_test
    LIBS parameter psi base device
    SOURCES test_transfer.cpp
            ../func_transfer.cpp
            ../base_matrix.cpp
            ../hcontainer.cpp
            ../atom_pair.cpp
            ../transfer.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp
            tmp_mocks.cpp)

  abacus_copy_test_files(parallel_hcontainer_tests.sh)
  find_program(BASH bash)
  add_test(
    NAME MODULE_LCAO_hcontainer_para_test
    COMMAND ${BASH} parallel_hcontainer_tests.sh
            "$<TARGET_FILE:MODULE_LCAO_hcontainer_transfer_test>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  addtest(
    TARGET MODULE_LCAO_hcontainer_output_test
    LIBS parameter base device
    SOURCES test_hcontainer_output.cpp
            tmp_mocks.cpp
            ../output_hcontainer.cpp
            ../base_matrix.cpp
            ../hcontainer.cpp
            ../atom_pair.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp)

  abacus_copy_test_files(support)

  addtest(
    TARGET MODULE_LCAO_hcontainer_add_value_test
    LIBS parameter ${math_libs} psi base device
    SOURCES test_add_value.cpp ../base_matrix.cpp ../hcontainer.cpp
            ../atom_pair.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp)

endif()
