if(ENABLE_LCAO)
  list(
    APPEND
    objects
    velocity_op.cpp
    pos_op_mat.cpp
    pos_op_basis.cpp
    pos_op_calc.cpp
    pos_op_writer.cpp
    pos_op_csr.cpp)

  add_library(momentum_op OBJECT ${objects})

  if(ENABLE_COVERAGE)
    add_coverage(momentum_op)
  endif()
endif()
