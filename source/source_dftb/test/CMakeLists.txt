AddTest(
  TARGET MODULE_DFTB_native
  LIBS parameter base device
  SOURCES test_eigensolver.cpp
    test_charge_mixing.cpp
    test_periodic_scc.cpp
    ../charge_mixing.cpp
    ../eigensolver.cpp
    ../periodic_scc.cpp
    ../sk_matrix.cpp
    ../sk_block.cpp
    ../skf_data.cpp
)
