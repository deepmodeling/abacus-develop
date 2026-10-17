if(NOT ENABLE_LIBRI OR NOT ENABLE_LCAO)
  return()
endif()

add_subdirectory(module_librpa_io)

add_library(rpa OBJECT
  rpa_lri.cpp
  rpa_lri_detail.cpp
  rpa_lri_band_io.cpp
  rpa_lri_coeff_v1.cpp
  rpa_lri_bz_io.cpp
  rpa_lri_mpi_layout.cpp
  rpa_lri_wfc_mpi.cpp
  rpa_lri_coulomb.cpp
  rpa_lri_coulomb_io.cpp
  rpa_lri_overlap.cpp
  rpa_lri_overlap_io.cpp
  rpa_lri_output.cpp
  rpa_lri_bz.cpp
  rpa_lri_sym.cpp
  rpa_lri_coeff_io.cpp
  rpa_lri_velocity.cpp
  rpa_lri_basis.cpp
  rpa_lri_eigenvector_v1.cpp
  rpa_lri_eigenvector.cpp
)
target_include_directories(rpa PRIVATE
  ${ABACUS_SOURCE_DIR}
  ${ABACUS_SOURCE_DIR}/source_lcao/module_ri
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/module_librpa_io
)
target_link_libraries(rpa PRIVATE ri librpa_io)

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(unittests)
endif()

if(ENABLE_COVERAGE)
  add_coverage(rpa)
endif()
