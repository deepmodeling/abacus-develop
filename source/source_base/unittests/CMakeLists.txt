abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES test_csr_file_reader.cpp test_csr_reader.cpp)

addtest(
  TARGET MODULE_BASE_filename
  LIBS parameter base device
  SOURCES test_filename.cpp)

addtest(
  TARGET MODULE_BASE_integral
  LIBS parameter base device
  SOURCES test_math_integral.cpp)

addtest(
  TARGET MODULE_BASE_net_force
  LIBS parameter base device
  SOURCES test_mathzone.cpp)
