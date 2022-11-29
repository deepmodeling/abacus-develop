list(APPEND operator_of_pw_srcs
    kedf_tf.cpp
    kedf_vw.cpp
    kedf_wt.cpp
)

add_library(
    operator_of_pw
    OBJECT
    ${operator_of_pw_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(operator_of_pw)
endif()
