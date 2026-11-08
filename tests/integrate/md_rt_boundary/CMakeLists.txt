# This contract intentionally spans MD and RT; keep it outside either unit target.
abacus_disable_feature_definitions(__LCAO)
abacus_add_local_feature_definitions(__NORMAL)

# Reuse the existing MD/cell support closure rather than duplicating its source list.
get_directory_property(
  md_test_sources DIRECTORY "${PROJECT_SOURCE_DIR}/source/source_md" DEFINITION
                            md_test_depend_files)

addtest(
  TARGET INTEGRATION_MD_RT_boundary
  LIBS parameter psi device container
  SOURCES
    test_md_rt_boundary.cpp
    ${PROJECT_SOURCE_DIR}/source/source_md/md_base.cpp
    ${PROJECT_SOURCE_DIR}/source/source_lcao/module_rt/boundary_fix.cpp
    ${PROJECT_SOURCE_DIR}/source/source_lcao/module_rt/td_folding.cpp
    ${PROJECT_SOURCE_DIR}/source/source_hamilt/module_hcontainer/base_matrix.cpp
    ${PROJECT_SOURCE_DIR}/source/source_hamilt/module_hcontainer/atom_pair.cpp
    ${PROJECT_SOURCE_DIR}/source/source_hamilt/module_hcontainer/hcontainer.cpp
    ${PROJECT_SOURCE_DIR}/source/source_basis/module_ao/parallel_orbitals.cpp
    ${PROJECT_SOURCE_DIR}/source/source_base/block_cyclic_distribution.cpp
    ${PROJECT_SOURCE_DIR}/source/source_base/parallel_common.cpp
    ${PROJECT_SOURCE_DIR}/source/source_cell/magnetism.cpp
    ${md_test_sources})
abacus_apply_build_options(INTEGRATION_MD_RT_boundary)

set_tests_properties(
  INTEGRATION_MD_RT_boundary
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;I_MPI_PIN=0"
             LABELS "integration;md;rt" TIMEOUT 60)
foreach(rank_count IN ITEMS 1 4)
  foreach(thread_count IN ITEMS 1 4)
    # AddTest already registers the singleton, single-thread case.
    if(rank_count EQUAL 1 AND thread_count EQUAL 1)
      continue()
    endif()
    set(test_name INTEGRATION_MD_RT_boundary_np${rank_count}_omp${thread_count})
    math(EXPR processor_count "${rank_count} * ${thread_count}")
    add_test(
      NAME ${test_name}
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count}
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:INTEGRATION_MD_RT_boundary>
        ${MPIEXEC_POSTFLAGS}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(
      ${test_name}
      PROPERTIES ENVIRONMENT
                 "OMP_NUM_THREADS=${thread_count};MKL_NUM_THREADS=1;I_MPI_PIN=0"
                 LABELS
                 "integration;md;rt"
                 PROCESSORS
                 ${processor_count}
                 TIMEOUT
                 60)
  endforeach()
endforeach()
