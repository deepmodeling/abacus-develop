if (ENABLE_LIBRI)
  AddTest(
    TARGET MODULE_LR_ri_benchmark_test
    LIBS psi base device container parameter
    SOURCES ri_benchmark_test.cpp
  )
endif()