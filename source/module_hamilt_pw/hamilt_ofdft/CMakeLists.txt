list(APPEND hamilt_ofdft_srcs
    kedf_tf.cpp
    kedf_vw.cpp
    kedf_wt.cpp
    kedf_lkt.cpp
    of_stress_pw.cpp
)

add_library(
    hamilt_ofdft
    OBJECT
    ${hamilt_ofdft_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_ofdft)
endif()
