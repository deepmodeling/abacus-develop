abacus_add_local_feature_definitions(__NORMAL)
if (USE_CUDA AND ENABLE_FLOAT_FFTW)
  AddTest(
    TARGET pw_test_gpu
    LIBS parameter base planewave device FFTW3::FFTW3_FLOAT
    SOURCES pw_test.cpp  pw_basis_C2R.cpp pw_basis_C2C.cpp  pw_basis_k_C2C.cpp 
  )
endif()

