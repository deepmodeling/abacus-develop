remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
if(ENABLE_LCAO)
  AddTest(
    TARGET EState_updaterhok_pw
    LIBS ${math_libs} planewave_serial base_serial psi
    SOURCES updaterhok_pw_test.cpp
            ../elecstate_pw.cpp ../elecstate.cpp
            ../../src_pw/charge.cpp ../../src_parallel/parallel_reduce.cpp
            ../../src_pw/structure_factor.cpp ../../src_pw/pw_complement.cpp
            ../../src_pw/klist.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_pw/occupy.cpp
            ../src/elecstate_multi_device.cpp
  )
  if(USE_CUDA)
    target_link_libraries(EState_updaterhok_pw cufft)
  endif()

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
              ../../src_lcao/local_orbital_charge.cpp ../../src_lcao/DM_gamma.cpp ../../src_lcao/DM_k.cpp
              ../../src_lcao/local_orbital_wfc.cpp
              ../../module_gint/grid_technique.cpp ../../module_gint/grid_meshball.cpp ../../module_gint/grid_bigcell.cpp
              ../../module_gint/grid_meshcell.cpp ../../module_gint/grid_meshk.cpp
              ../../module_gint/gint_tools.cpp ../../module_gint/gint_k_pvpr.cpp ../../module_gint/gint_tau.cpp
              ../../src_lcao/dm_2d.cpp ../../src_lcao/LCAO_matrix.cpp
              ../../src_lcao/record_adj.cpp  ../../src_lcao/LCAO_nnr.cpp
              ../../module_gint/gint.cpp ../../module_gint/gint_gamma.cpp
              ../../module_gint/gint_fvl.cpp ../../module_gint/gint_rho.cpp
              ../../module_gint/gint_vl.cpp
              ../../src_pw/charge.cpp
              ../../src_pdiag/pdiag_common.cpp
              ../../src_io/output.cpp ../../src_pw/soc.cpp ../../src_io/read_rho.cpp
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
add_library(
  base_serial
  OBJECT
  ../../module_base/complexarray.cpp
  ../../module_base/complexmatrix.cpp
  ../../module_base/element_basis_index.cpp
  ../../module_base/export.cpp
  ../../module_base/integral.cpp
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
  ../../module_base/mathzone.cpp
  ../../module_base/mathzone_add1.cpp
  ../../module_base/matrix.cpp
  ../../module_base/matrix3.cpp
  ../../module_base/memory.cpp
  ../../module_base/mymath1.cpp
  ../../module_base/mymath3.cpp
  ../../module_base/polint.cpp
  ../../module_base/realarray.cpp
  ../../module_base/sph_bessel.cpp
  ../../module_base/sph_bessel_recursive-d1.cpp
  ../../module_base/sph_bessel_recursive-d2.cpp
  ../../module_base/timer.cpp
  ../../module_base/tool_check.cpp
  ../../module_base/tool_quit.cpp
  ../../module_base/tool_title.cpp
  ../../module_base/ylm.cpp
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
  ../../module_pw/src/pw_multi_device.cpp
)

if (USE_CUDA)
  list(APPEND planewave_serial_srcs ../../module_pw/src/cuda/pw_multi_device.cu)
endif()

add_library(
  planewave_serial
  OBJECT
  ${planewave_serial_srcs}
)

if (USE_CUDA)
AddTest(
  TARGET Elecstate_UTs
  LIBS ${math_libs} psi base
  SOURCES elecstate_multi_device_test.cpp
          ../src/elecstate_multi_device.cpp
          ../src/cuda/elecstate_multi_device.cu
          ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/parallel_global.cpp 
)
else()
AddTest(
  TARGET Elecstate_UTs
  LIBS ${math_libs} psi base
  SOURCES elecstate_multi_device_test.cpp
          ../src/elecstate_multi_device.cpp
          ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/parallel_global.cpp 
)
endif()

if(ENABLE_COVERAGE)
  add_coverage(base_serial)
  add_coverage(planewave_serial)
endif()
