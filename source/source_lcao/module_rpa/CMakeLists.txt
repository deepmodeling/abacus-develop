if(NOT ENABLE_LIBRI)
  return()
endif()

add_library(rpa OBJECT
  rpa_abfs_preorthogonalization.cpp
)
target_include_directories(rpa
  PRIVATE
    ${ABACUS_SOURCE_DIR}
    ${ABACUS_SOURCE_DIR}/source_lcao/module_ri
    ${ABACUS_SOURCE_DIR}/source_lcao/module_rpa
    ${ABACUS_SOURCE_DIR}/source_lcao/module_rpa/module_librpa_io
)
target_link_libraries(rpa PRIVATE ri)

add_subdirectory(module_delta_st)
add_subdirectory(module_librpa_io)

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(test)
endif()

if(ENABLE_COVERAGE)
  add_coverage(rpa)
endif()
