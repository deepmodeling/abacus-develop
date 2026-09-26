# Restricted rVV10 caller: exercise the same compiled objects as the executable.
if(NOT USE_CUDA AND NOT USE_ROCM)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME RVV10_DRIVER_OUTPUT_VALIDATION
    COMMAND ${Python3_EXECUTABLE} -m unittest discover
      -s ${CMAKE_CURRENT_SOURCE_DIR} -p test_regression_driver.py)
  if(ENABLE_MPI AND ENABLE_LIBXC)
    set(rvv10_modes scf spin-scf mpi-scf mpi-spin-scf mpi-spin-reject mpi-pbe)
  elseif(ENABLE_LIBXC)
    set(rvv10_modes scf spin-scf reject uspp)
  else()
    set(rvv10_modes no-libxc)
  endif()
  foreach(mode IN LISTS rvv10_modes)
    set(rvv10_launcher)
    set(rvv10_ranks 1)
    if(ENABLE_MPI)
      if(mode MATCHES "^mpi-")
        set(rvv10_ranks 2)
      endif()
      set(rvv10_launcher --launch ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
        ${rvv10_ranks} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${ABACUS_BIN_NAME}> ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(NAME RVV10_${mode}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_scf_regression.py
        --executable $<TARGET_FILE:${ABACUS_BIN_NAME}>
        --pseudo-dir ${PROJECT_SOURCE_DIR}/tests/PP_ORB
        --artifacts-dir ${CMAKE_CURRENT_BINARY_DIR}/artifacts --mode ${mode} ${rvv10_launcher})
    set_tests_properties(RVV10_${mode} PROPERTIES
      TIMEOUT 1700 PROCESSORS ${rvv10_ranks}
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
  endforeach()
  if(ENABLE_LIBXC AND NOT ENABLE_MPI)
    get_target_property(rvv10_native_libraries ${ABACUS_BIN_NAME} LINK_LIBRARIES)
    AddTest(
      TARGET MODULE_ESTATE_RVV10_POTENTIAL
      SOURCES test_pot_rvv10.cpp
      LIBS ${rvv10_native_libraries}
    )
    abacus_apply_build_options(MODULE_ESTATE_RVV10_POTENTIAL)
    set_tests_properties(MODULE_ESTATE_RVV10_POTENTIAL PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
  endif()
endif()
