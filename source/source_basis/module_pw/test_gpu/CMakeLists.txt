abacus_add_local_feature_definitions(__NORMAL)
if (USE_CUDA AND ENABLE_FLOAT_FFTW)
  AddTest(
    TARGET pw_test_gpu
    LIBS parameter base planewave device FFTW3::FFTW3_FLOAT
    SOURCES pw_test.cpp  pw_basis_c2r.cpp pw_basis_c2c.cpp  pw_basis_k_c2c.cpp
  )
endif()

