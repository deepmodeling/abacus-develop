abacus_disable_feature_definitions(__LIBXC)
addtest(
  TARGET MODULE_LR_lr_util_phys_test
  LIBS parameter base device container planewave #for FFT
  SOURCES lr_util_physics_test.cpp ../lr_util.cpp)

addtest(
  TARGET MODULE_LR_lr_util_algo_test
  LIBS parameter base device psi container planewave #for FFT
  SOURCES lr_util_algo_test.cpp ../lr_util.cpp)

addtest(
  TARGET MODULE_LR_lr_util_pp
  LIBS parameter
       base
       device
       container
       planewave
       cell_info
       symmetry
       psi
  SOURCES test_lr_util.cpp ../lr_util.cpp ../../../../source_cell/magnetism.cpp)

addtest(
  TARGET MODULE_LR_spectrum_mo_output_test
  LIBS parameter base device psi container planewave
  SOURCES test_spectrum_mo_output.cpp ../spectrum_out.cpp)
