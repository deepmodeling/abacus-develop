AddTest(
  TARGET MODULE_DFTB_native_eigensolver
  LIBS parameter base device
  SOURCES native_eigensolver_test.cpp
    ../eigensolver.cpp
    ../sk_matrix.cpp
    ../sk_block.cpp
    ../skf_data.cpp
)
