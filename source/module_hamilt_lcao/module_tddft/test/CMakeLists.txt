remove_definitions(-D __MPI)

AddTest(
  TARGET tddft_test
  LIBS ${math_libs} base device
  SOURCES tddft_test.cpp upsi_test1.cpp upsi_test2.cpp upsi_test3.cpp propagator_test1.cpp propagator_test2.cpp propagator_test3.cpp middle_hamilt_test.cpp norm_psi_test.cpp bandenergy_test.cpp ../upsi.cpp ../propagator.cpp ../middle_hamilt.cpp ../norm_psi.cpp ../bandenergy.cpp 
)

add_test(NAME tddft_test_np1
      COMMAND mpirun -np 1 ./tddft_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

