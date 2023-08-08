remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)
remove_definitions(-D_OPENMP)

if (USE_ABACUS_LIBM)
list (APPEND LIBM_SRC_serial
  ../../module_base/libm/branred.cpp
  ../../module_base/libm/cexp.cpp
  ../../module_base/libm/exp.cpp
  ../../module_base/libm/sincos.cpp
)
endif()

add_library(
    base_serial
    OBJECT
    ../../module_base/complexarray.cpp
    ../../module_base/complexmatrix.cpp
    ../../module_base/element_basis_index.cpp
    ../../module_base/export.cpp
    ../../module_base/inverse_matrix.cpp
    ../../module_base/gather_math_lib_info.cpp
    ../../module_base/global_file.cpp
    ../../module_base/global_function.cpp
    ../../module_base/global_function_ddotreal.cpp
    ../../module_base/global_variable.cpp
    ../../module_base/intarray.cpp
    ../../module_base/math_integral.cpp
    ../../module_base/math_polyint.cpp
    ../../module_base/math_sphbes.cpp
    ../../module_base/math_ylmreal.cpp
    ../../module_base/math_bspline.cpp
    ../../module_base/math_chebyshev.cpp
    ../../module_base/mathzone_add1.cpp
    ../../module_base/matrix.cpp
    ../../module_base/matrix3.cpp
    ../../module_base/memory.cpp
    ../../module_base/mymath.cpp
    ../../module_base/opt_CG.cpp
    ../../module_base/opt_DCsrch.cpp
    ../../module_base/realarray.cpp
    ../../module_base/sph_bessel_recursive-d1.cpp
    ../../module_base/sph_bessel_recursive-d2.cpp
    ../../module_base/timer.cpp
    ../../module_base/tool_check.cpp
    ../../module_base/tool_quit.cpp
    ../../module_base/tool_title.cpp
    ../../module_base/ylm.cpp
    ../../module_base/abfs-vector3_order.cpp
    ../../module_base/parallel_common.cpp
    ../../module_base/parallel_global.cpp
    ../../module_base/parallel_reduce.cpp
    ../../module_base/spherical_bessel_transformer.cpp
    ../../module_base/cubic_spline.cpp
    ${LIBM_SRC_serial}
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET io_rho_io
  LIBS ${math_libs} base device cell_info
  SOURCES rho_io_test.cpp ../read_cube.cpp ../write_cube.cpp ../read_rho.cpp ../write_rho.cpp ../output.cpp
)

AddTest(
  TARGET io_dm_io
  LIBS ${math_libs} base_serial device
  SOURCES dm_io_test.cpp ../read_dm.cpp ../write_dm.cpp ../output.cpp
  ../../module_cell/unitcell.cpp
  ../../module_cell/read_atoms.cpp
  ../../module_cell/read_cell_pseudopots.cpp
  ../../module_cell/atom_spec.cpp
  ../../module_cell/atom_pseudo.cpp
  ../../module_cell/pseudo_nc.cpp
  ../../module_cell/read_pp.cpp
  ../../module_cell/read_pp_upf201.cpp
  ../../module_cell/read_pp_upf100.cpp
  ../../module_cell/read_pp_vwr.cpp
  ../../module_cell/read_pp_blps.cpp
)

AddTest(
  TARGET io_nscf_band
  LIBS ${math_libs} base device
  SOURCES nscf_band_test.cpp ../nscf_band.cpp
)

AddTest(
  TARGET io_read_wfc_nao_test
  LIBS ${math_libs} base device
  SOURCES read_wfc_nao_test.cpp ../read_wfc_nao.cpp ../../module_psi/psi.cpp
)
