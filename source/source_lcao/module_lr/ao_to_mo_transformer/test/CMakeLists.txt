abacus_disable_feature_definitions(__LIBXC)
AddTest(
  TARGET MODULE_LR_ao_to_mo_test
  LIBS parameter base container device psi
  SOURCES ao_to_mo_test.cpp ../../utils/lr_util.cpp ../ao_to_mo_parallel.cpp ../ao_to_mo_serial.cpp
)

AddTest(
  TARGET MODULE_LR_CVCX_test
  LIBS base parameter ${math_libs} container device psi
  SOURCES test_cvcx.cpp ../../utils/lr_util.cpp ../cvcx_par.cpp ../cvcx_serial.cpp
)