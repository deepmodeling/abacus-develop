if(ENABLE_LCAO)
  list(
    APPEND
    objects
    velocity_op.cpp
    pos_op_mat.cpp
    pos_op_basis.cpp
    pos_op_calc.cpp
    pos_op_writer.cpp
    pos_op_csr.cpp
    mom_op_mat.cpp
    angmom_op_mat.cpp)

  add_library(op OBJECT ${objects})

  if(ENABLE_COVERAGE)
    add_coverage(op)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(unittests)
  endif()
endif()
