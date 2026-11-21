abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES test_csr_file_reader.cpp test_csr_reader.cpp)

AddTest(
  TARGET MODULE_BASE_filename
  LIBS parameter base device
  SOURCES test_filename.cpp
)
