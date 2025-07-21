remove_definitions(-D__MPI)

AddTest(
  TARGET MODULE_BASE_GRID_test_delley
  SOURCES test_delley.cpp
  ../delley.cpp
  ../../ylm.cpp
)

AddTest(
  TARGET MODULE_BASE_GRID_test_radial
  SOURCES test_radial.cpp
  ../radial.cpp
)

AddTest(
  TARGET MODULE_BASE_GRID_test_partition
  SOURCES test_partition.cpp
  ../partition.cpp
  ../radial.cpp
  ../delley.cpp
)

AddTest(
  TARGET MODULE_BASE_GRID_test_batch
  SOURCES test_batch.cpp
  ../batch.cpp
  ../../blas_connector_base.cpp
  ../../blas_connector_l1.cpp
  ../../blas_connector_l2.cpp
  ../../blas_connector_l3.cpp
  ../../lapack_connector.cpp
)
