# These tests exercise CPU code and need no GPU kernels.
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

set(ofdft_test_sources test_kedf_tf.cpp ../kedf_tf.cpp)
set(ofdft_test_libs parameter base device)
if(ENABLE_MPI)
  # The linked base library performs MPI reductions even for a single rank.
  list(APPEND ofdft_test_sources ../../../source_base/test/mpi_test_main.cpp)
  list(APPEND ofdft_test_libs MPI::MPI_CXX)
endif()

AddTest(
  TARGET MODULE_PW_OfDFT
  LIBS ${ofdft_test_libs}
  SOURCES ${ofdft_test_sources}
)
