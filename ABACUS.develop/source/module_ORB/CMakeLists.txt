add_library(
    orb
    OBJECT
    ORB_control.cpp
    ORB_read.cpp
    ORB_atomic.cpp
    ORB_atomic_lm.cpp
    ORB_nonlocal.cpp
    ORB_nonlocal_lm.cpp
    ORB_gaunt_table.cpp
    ORB_table_beta.cpp
    ORB_table_phi.cpp
    ORB_table_alpha.cpp
    ORB_gen_tables.cpp
)
target_compile_definitions(orb PRIVATE __NORMAL)
target_compile_options(
    orb
    PRIVATE
    -Ofast -std=c++11 -simd -xHost -m64 -qopenmp -Werror -Wall -pedantic -g
)
target_link_libraries(orb PRIVATE global)
