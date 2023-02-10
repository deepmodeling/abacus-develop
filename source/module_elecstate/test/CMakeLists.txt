remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)

if(ENABLE_LCAO)
  AddTest(
    TARGET EState_updaterhok_pw
    LIBS ${math_libs} planewave_serial base_serial psi device
    SOURCES updaterhok_pw_test.cpp
            ../elecstate_pw.cpp ../elecstate.cpp
            ../module_charge/charge.cpp ../../src_parallel/parallel_reduce.cpp
            ../../module_hamilt_pw/hamilt_pwdft/structure_factor.cpp
            ../../module_cell/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../occupy.cpp
            ../kernels/elecstate_op.cpp ../../module_io/rwstream.cpp ../../module_io/write_wfc_pw.cpp
  )

  install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  list(APPEND test_math_libs ${math_libs})
  if(USE_ELPA)
    list(APPEND test_math_libs ELPA::ELPA)
  endif()
  #[[AddTest(
    TARGET EState_psiToRho_lcao
    LIBS ${test_math_libs} base orb cell neighbor planewave psi
    SOURCES elecstate_lcao_test.cpp ../elecstate_lcao.cpp
              ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp
              ../../module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.cpp ../../module_hamilt_lcao/hamilt_lcaodft/DM_gamma.cpp ../../module_hamilt_lcao/hamilt_lcaodft/DM_k.cpp
              ../../module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.cpp
              ../../module_gint/grid_technique.cpp ../../module_gint/grid_meshball.cpp ../../module_gint/grid_bigcell.cpp
              ../../module_gint/grid_meshcell.cpp ../../module_gint/grid_meshk.cpp
              ../../module_gint/gint_tools.cpp ../../module_gint/gint_k_pvpr.cpp ../../module_gint/gint_tau.cpp
              ../../module_hamilt_lcao/hamilt_lcaodft/dm_2d.cpp ../../module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.cpp
              ../../module_hamilt_lcao/hamilt_lcaodft/record_adj.cpp  ../../module_hamilt_lcao/hamilt_lcaodft/LCAO_nnr.cpp
              ../../module_gint/gint.cpp ../../module_gint/gint_gamma.cpp
              ../../module_gint/gint_fvl.cpp ../../module_gint/gint_rho.cpp
              ../../module_gint/gint_vl.cpp
              ../../module_hamilt_pw/hamilt_pwdft/charge.cpp
              ../../module_io/output.cpp ../../module_hamilt_pw/hamilt_pwdft/soc.cpp ../../module_io/read_rho.cpp
  )
  if(USE_CUDA)
    target_link_libraries(EState_psiToRho_lcao cufft)
  endif()
  target_compile_definitions(EState_psiToRho_lcao PRIVATE __MPI)
  install(FILES elecstate_lcao_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  find_program(BASH bash)
  add_test(NAME EState_psiToRho_lcao_parallel
        COMMAND ${BASH} elecstate_lcao_parallel_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )]]
endif()
if (USE_ABACUS_LIBM)
list (APPEND LIBM_SRC
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
  ../../module_base/realarray.cpp
  ../../module_base/sph_bessel_recursive-d1.cpp
  ../../module_base/sph_bessel_recursive-d2.cpp
  ../../module_base/timer.cpp
  ../../module_base/tool_check.cpp
  ../../module_base/tool_quit.cpp
  ../../module_base/tool_title.cpp
  ../../module_base/ylm.cpp
  ${LIBM_SRC}
)

list(APPEND planewave_serial_srcs
  ../../module_pw/fft.cpp
  ../../module_pw/pw_basis.cpp
  ../../module_pw/pw_basis_k.cpp
  ../../module_pw/pw_distributeg.cpp
  ../../module_pw/pw_distributeg_method1.cpp
  ../../module_pw/pw_distributeg_method2.cpp
  ../../module_pw/pw_distributer.cpp
  ../../module_pw/pw_init.cpp
  ../../module_pw/pw_transform.cpp
  ../../module_pw/pw_transform_k.cpp
)

add_library(
  planewave_serial
  OBJECT
  ${planewave_serial_srcs}
)

AddTest(
  TARGET Elecstate_UTs
  LIBS ${math_libs} psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
          ../../src_parallel/parallel_reduce.cpp
          ../../src_parallel/parallel_global.cpp 
)

if(ENABLE_COVERAGE)
  add_coverage(base_serial)
  add_coverage(planewave_serial)
endif()
