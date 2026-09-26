add_library(abacus_greenx_minimax STATIC
  GX-common/src/kinds.f90
  GX-common/src/constants.f90
  GX-common/src/unit_conversion.f90
  GX-common/src/error_handling.f90
  GX-TimeFrequency/src/minimax_utils.F90
  GX-TimeFrequency/src/minimax_tau.F90
  GX-TimeFrequency/src/minimax_omega.F90
  GX-TimeFrequency/src/minimax_grids.F90
  GX-TimeFrequency/api/api_utilities.f90
  GX-TimeFrequency/api/gx_minimax.f90
  gx_minimax_wrp.f90
)

set_target_properties(abacus_greenx_minimax PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# The bundled GreenX tables intentionally contain long free-form Fortran
# literals.  GNU Fortran otherwise truncates them at column 132 (and the CI
# treats that warning as an error), while Intel Fortran accepts the sources as
# written.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(abacus_greenx_minimax PRIVATE -ffree-line-length-none)
endif()

target_include_directories(abacus_greenx_minimax PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(abacus_greenx_minimax PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/GX-common/src
  ${CMAKE_CURRENT_SOURCE_DIR}/GX-TimeFrequency/src
  ${CMAKE_CURRENT_SOURCE_DIR}/GX-TimeFrequency/api
)

target_link_libraries(abacus_greenx_minimax PUBLIC ${math_libs})
