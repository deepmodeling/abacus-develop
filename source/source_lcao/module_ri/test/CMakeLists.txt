abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
if(ENABLE_LCAO)
  AddTest(
    TARGET MODULE_RI_exx_lr_ws
    LIBS base parameter device ${math_libs}
    SOURCES test_exx_lr_ws.cpp ../exx_lr_ws.cpp ../../../source_basis/module_ao/orb_atomic_lm.cpp
  )
endif()
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES dm_mixing_test.cpp ../mix_dmk_2d.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES ri_cv_io_test.cpp
)
AddTest(
  TARGET MODULE_RI_abfs_vector3_order
  LIBS parameter
  SOURCES abfs_vec3_order_test.cpp
)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
