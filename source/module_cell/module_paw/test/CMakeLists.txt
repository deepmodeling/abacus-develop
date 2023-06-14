remove_definitions(-D__MPI)

install(FILES H.LDA_PW-JTH.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET test_paw
  LIBS ${math_libs} base device
  SOURCES test_read_paw.cpp ../paw_element.cpp
)