find_package(Python3 QUIET COMPONENTS Interpreter)

if(Python3_Interpreter_FOUND)
  add_test(NAME RVV10_LIBVDWXC_CONTRACT
    COMMAND ${Python3_EXECUTABLE} -m unittest discover
      -s ${CMAKE_CURRENT_SOURCE_DIR} -p test_libvdwxc_contract.py)

  set_tests_properties(RVV10_LIBVDWXC_CONTRACT PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1")
endif()
